# MIT License

# Copyright (c) 2024 Samuel Waldner s.waldner@unibas.ch

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
import pandas as pd
import skimage
import tifffile as tif
import os
import shutil
import openpyxl
from glob import glob
import matplotlib.pyplot as plt
import configparser

#set default values
version = 1.00
data_dir = os.path.dirname(__file__)
nucleus_area_cutoff = 7
overlap_area_cutoff = 7
z_dimension_scale_factor = 0.1035
verbose = False
create_histogram = False
count_green_puncta = False

### Input Parameters ###
config = configparser.ConfigParser()
config.read('settings.ini')
nucleus_area_cutoff = float(config['Settings']['Nucleus_cutoff_sidelength'])            #sidelength of the square (for 2D) or the cube (in 3D) that defines the minimal acceptable area (in pixels) of a nucleus. Nuclei smaller than this are removed before analysis
overlap_area_cutoff = float(config['Settings']['Overlap_cutoff_sidelength'])            #sidelength of the square (for 2D) or the cube (in 3D) that defines the maximum acceptable sized of colocalized areas (in pixels) of blue and green. Spots larger than this are removed before analysis
z_dimension_scale_factor = float(config['Settings']['z_dimension_scale_factor'])        #how much smaller the scale of the x/y-axis is (for z-stacks) than of the z-axis -> is used to adequately calculate the overlap/nucleus cutoff values for 3D
if config['Settings']['Save_intermediate_images'] == 'True':                            #if set to true, creates subdirectory and saves overlap .tif and labelled cell .tif images. Useful for validating results/bug fixing/visualization
    verbose = True
elif config['Settings']['Save_intermediate_images'] == 'False':
    verbose = False
if config['Settings']['Create_histograms'] == 'True':                                   #if set to true, creates subdirectory and saves histograms of events per cell data per timepoint as png images. Useful for validation of the quantification as well as general science
    create_histogram = True
elif config['Settings']['Create_histograms'] == 'False':
    create_histogram = False
if config['Settings']['Count_green_puncta'] == 'True':                                  #if set to true, counts green signal puncta in addition to overlap (green-blue) puncta.
    count_green_puncta = True
elif config['Settings']['Count_green_puncta'] == 'False':
    count_green_puncta = False                   
########################

#loop through subdirectories

def process_directory(directory):
    parent_folder_name = os.path.split(directory)[1]
    # Check if '*.xlsx' file already exists in the directory
    xlsx_path = os.path.join(directory, f'{parent_folder_name}_results.xlsx')
    existing_xlsx_files = glob(xlsx_path)

    if existing_xlsx_files:
        print(f"Excel file already exists in {directory}")
    else:
        # Look for '*.tiff' files containing 'segmented' in their filenames
        tiff_files = glob(os.path.join(directory, '*segmented*.tiff'))
        # Check if there are exactly two '*segmented*.tiff' files
        if len(tiff_files) == 4:
            # Copy the template '*.xlsx' file into the directory
            copy_template_xlsx(directory,parent_folder_name)

            # Process the contents of the '*.tiff' files and write to '*.xlsx'
            if process_tiff_files(tiff_files, directory, parent_folder_name) == True:
                print(f"Results.xlsx file created in {directory}")
        else:
            print(f"Error: Four '*segmented*.tiff' files not found in {directory}")

def copy_template_xlsx(directory,parent_folder_name):
    # Specify the path to the template Excel file
    template_path = os.path.join(data_dir,'template_per_cell_data.xlsx')

    # Specify the destination path for the new Excel file
    destination_path = os.path.join(directory, f'{parent_folder_name}_results.xlsx')

    # Copy the template Excel file to the destination path
    shutil.copy(template_path, destination_path)

def process_tiff_files(tiff_files, directory,parent_folder_name):
    #if verbose mode is active create directory to save intermediary images into:
    if not os.path.exists(os.path.join(directory, 'labeleled_cells')) and verbose:
        try:
            os.makedirs(os.path.join(directory, 'labeleled_cells'))
            os.makedirs(os.path.join(directory, 'labeled_puncta'))
            os.makedirs(os.path.join(directory, 'labeled_colocalized_spots'))
        except:
            print('mkdir labeled_cells failed')
            return(False)

    if not os.path.exists(os.path.join(directory, 'histograms')) and create_histogram:
        try:
            os.makedirs(os.path.join(directory, 'histograms'))
        except:
            print('mkdir histogram failed')
            return(False)        

    #read tiff files
    blue_loaded = green_loaded = nuclei_loaded = cells_loaded = False
    for tiff_file in tiff_files:
        if 'blue' in tiff_file or 'blu3' in tiff_file:
            image_blue_segmented = tif.imread(tiff_file)

            image_blue_path_head,image_blue_path_tail = os.path.split(tiff_file)
            image_blue_path_tail = image_blue_path_tail.replace('_segmented','')
            image_blue_path_tail = image_blue_path_tail.replace('.tiff','.tif')
            image_blue_path = os.path.join(image_blue_path_head,image_blue_path_tail)
            if os.path.isfile(image_blue_path):
                image_blue = tif.imread(image_blue_path)
            else:
                print('Unable to find labelled DNA image file unsegmented version - stopping script')
                return(False)
                #quit()
            blue_loaded = True
        elif ('green' in tiff_file or 'gr33n' in tiff_file) and ('cells' not in tiff_file):
            image_green_segmented = tif.imread(tiff_file)
            green_loaded = True
        elif 'nuclei' in tiff_file or 'nucl3i' in tiff_file:
            image_nuclei_segmented = tif.imread(tiff_file)
            nuclei_loaded = True
        elif 'cells' in tiff_file:
            image_cells_segmented = tif.imread(tiff_file)
            cells_loaded = True
        else:
            print('Unexpected filename encountered, stopping script.')
            return(False)
            #quit()
        
    if not (blue_loaded and green_loaded and nuclei_loaded and cells_loaded):
        print(f'Not all files loaded, stopping script. Blue: {blue_loaded}, Green: {green_loaded}, Nuclei: {nuclei_loaded}, Cells: {cells_loaded}')
        return(False)
        #quit()

    image_blue_segmented[image_blue_segmented == 1] = 0
    image_blue_segmented[image_blue_segmented == 2] = 1
    image_green_segmented[image_green_segmented == 1] = 0
    image_green_segmented[image_green_segmented == 2] = 1
    image_cells_segmented[image_cells_segmented == 1] = 0
    image_cells_segmented[image_cells_segmented == 2] = 1
    image_nuclei_segmented[image_nuclei_segmented == 1] = 0

    #multiply segmentation green & segmentation blue & segmentation cells
    image_overlap = np.multiply(image_green_segmented, image_blue_segmented)
    image_overlap = np.multiply(image_overlap, image_cells_segmented)

    image_DNA_in_cells = np.copy(image_blue)
    image_DNA_in_cells[image_cells_segmented == 0] = 0
    image_DNA_outside_cells = np.copy(image_blue)
    image_DNA_outside_cells[image_cells_segmented != 0] = 0

    image_dimensions = image_overlap.ndim
    if image_dimensions == 3:
        t_max,y_max,x_max = image_overlap.shape
    elif image_dimensions == 4:
        t_max,z_max,y_max,x_max = image_overlap.shape
    time_iterator = 0
    overlap_df = pd.DataFrame(columns=['Time', 'Overlap_label', 'Cell_label'])
    green_puncta_df = pd.DataFrame(columns=['Time', 'Overlap_label', 'Cell_label'])
    cell_df = pd.DataFrame(columns=['Time', 'Cell_label', 'Cell_area', 'DNA_area', 'Overlap_area', 'Overlap_area_percent_of_DNA_area', 'DNA_mean_signal_intensity'])
    dna_count_df = pd.DataFrame(columns=['Time','DNA_in_cells','DNA_outside_cells','area_inside_cells','area_outside_cells','DNA_in_cells_per_area','DNA_outside_cells_per_area'])
    
    while time_iterator < t_max:
        print(f'working on timepoint {time_iterator:03} out of {t_max:03} - [ ] Cleaning nucleus segmentation image', end='\r')

        #clean nucleus images: keep only large enough blobs
        if image_dimensions == 3:
            image_cells_segmented_current = np.copy(image_cells_segmented[time_iterator,:,:])
            image_nuclei_segmented_current = np.copy(image_nuclei_segmented[time_iterator,:,:])
            image_DNA_in_cells_current = np.copy(image_DNA_in_cells[time_iterator,:,:])
            image_DNA_outside_cells_current = np.copy(image_DNA_outside_cells[time_iterator,:,:])
            image_overlap_current = np.copy(image_overlap[time_iterator,:,:])
            image_blue_segmented_current = np.copy(image_blue_segmented[time_iterator,:,:])
            image_green_segmented_current = np.copy(image_green_segmented[time_iterator,:,:])
            overlap_area_cutoff_adjusted = overlap_area_cutoff * overlap_area_cutoff
            nucleus_area_cutoff_adjusted = nucleus_area_cutoff * nucleus_area_cutoff
        elif image_dimensions == 4:
            image_cells_segmented_current = np.copy(image_cells_segmented[time_iterator,:,:,:])
            image_nuclei_segmented_current = np.copy(image_nuclei_segmented[time_iterator,:,:,:])
            image_DNA_in_cells_current = np.copy(image_DNA_in_cells[time_iterator,:,:,:])
            image_DNA_outside_cells_current = np.copy(image_DNA_outside_cells[time_iterator,:,:,:])
            image_overlap_current = np.copy(image_overlap[time_iterator,:,:,:])
            image_blue_segmented_current = np.copy(image_blue_segmented[time_iterator,:,:,:])
            image_green_segmented_current = np.copy(image_green_segmented[time_iterator,:,:,:])
            overlap_area_cutoff_adjusted = overlap_area_cutoff * overlap_area_cutoff * (overlap_area_cutoff * z_dimension_scale_factor)
            nucleus_area_cutoff_adjusted = nucleus_area_cutoff * nucleus_area_cutoff * (nucleus_area_cutoff * z_dimension_scale_factor)
            
        #label individual nuclei
        image_nuclei_segmented_labels = np.copy(image_nuclei_segmented_current)
        image_nuclei_segmented_labels = skimage.measure.label(image_nuclei_segmented_labels, background=0)
        
        #prepare list of nuclei to be removed from the count due to their small size
        nuclei_list_to_clean = []

        #count area of each nucleus, fill in list to be cleaned
        moments = skimage.measure.regionprops_table(image_nuclei_segmented_labels, properties=('label','area'))
        current_cell_count = moments['label'][-1]
        for current_label in moments['label']:
            if moments['area'][current_label-1] < nucleus_area_cutoff_adjusted:
                nuclei_list_to_clean.append(current_label)
                current_area = moments['area'][current_label-1]
                current_cell_count -= 1

        #do actual cleaning on the original nucleus segmentation image
        for cleaning_label in nuclei_list_to_clean:
            image_nuclei_segmented_current[image_nuclei_segmented_labels == cleaning_label] = 0

        #label that image
        image_nuclei_segmented_current = skimage.measure.label(image_nuclei_segmented_current, background=0)

        print(f'working on timepoint {time_iterator:03} out of {t_max:03} - [x] Cleaning nucleus segmentation image - [ ] Segmenting individual cells',end='\r')
        #perform watershed segmentation using labelled nuclei as troughs, resulting in separated cells
        image_cells_segmented_current = skimage.segmentation.watershed(image_cells_segmented_current, markers=image_nuclei_segmented_current, mask=(image_cells_segmented_current > 0))

        print(f'working on timepoint {time_iterator:03} out of {t_max:03} - Quantifying DNA outside the cells      ',end='\r')

        #count DNA (by summing intesity) outside all cells. Normalize by area
        current_DNA_count_outside_cells = np.sum(image_DNA_outside_cells_current)
        current_area_outside_cells = (image_cells_segmented_current == 0).sum()
        current_DNA_count_outside_cells_by_area = current_DNA_count_outside_cells / current_area_outside_cells


        #remove cells that touch the image edge:
        if image_dimensions == 3:
            edge_1 = image_cells_segmented_current[0,:]
            edge_2 = image_cells_segmented_current[-1,:]
            edge_3 = image_cells_segmented_current[:,0]
            edge_4 = image_cells_segmented_current[:,-1]
            edge_total = np.concatenate((edge_1,edge_2,edge_3,edge_4))
        elif image_dimensions == 4: #in z-stack images ignore top and bottom edge (as most/all cells in the current dataset touch those edges) Uncommenting 'edge_1' and 'edge_2' lines and swapping over the np.concatenate command removes that restriction
            #edge_1 = image_cells_segmented_current[0,:,:]
            #edge_2 = image_cells_segmented_current[-1,:,:]
            edge_3 = image_cells_segmented_current[:,0,:]
            edge_4 = image_cells_segmented_current[:,-1,:]
            edge_5 = image_cells_segmented_current[:,:,0]
            edge_6 = image_cells_segmented_current[:,:,-1]
            edge_total = np.concatenate((edge_3,edge_4,edge_5,edge_6))#np.concatenate((edge_1,edge_2,edge_3,edge_4,edge_5,edge_6))
        edge_total_unique = np.unique(edge_total)

        print(f'working on timepoint {time_iterator:03} out of {t_max:03} - [x] Cleaning nucleus segmentation image - [x] Segmenting individual cells - [ ] Removing {len(edge_total_unique)} cells that touch the edge',end='\r')

        for element in edge_total_unique:
            image_cells_segmented_current[image_cells_segmented_current == element] = 0

        print(f'working on timepoint {time_iterator:03} out of {t_max:03} - [x] Cleaning nucleus segmentation image - [x] Segmenting individual cells - [x] Removing {len(edge_total_unique)} cells that touch the edge - [ ] Quantifying DNA inside the cells',end='\r')

        #count DNA (by summing intesity) inside all cells. Normalize by area
        image_DNA_in_cells_current[image_cells_segmented_current == 0] = 0
        current_DNA_count_in_cells = np.sum(image_DNA_in_cells_current)
        current_area_in_cells = (image_cells_segmented_current > 0).sum()
        current_DNA_count_in_cells_by_area = current_DNA_count_in_cells / current_area_in_cells
        dna_count_df.loc[len(dna_count_df)] = [int(time_iterator), current_DNA_count_in_cells, current_DNA_count_outside_cells, current_area_in_cells, current_area_outside_cells, current_DNA_count_in_cells_by_area, current_DNA_count_outside_cells_by_area]
        
        print(f'working on timepoint {time_iterator:03} out of {t_max:03} - [x] Cleaning nucleus segmentation image - [x] Segmenting individual cells - [x] Removing {len(edge_total_unique)} cells that touch the edge - [x] Quantifying DNA inside the cells - [ ] Characterizing individual cells',end='\r')
        
        image_cells_segmented_current_unique, image_cells_segmented_current_counts = np.unique(image_cells_segmented_current, return_counts=True)
        
        max_iterator = len(image_cells_segmented_current_unique)
        for iterator in range(max_iterator):
            if iterator > 0:
                current_label = image_cells_segmented_current_unique[iterator]
                current_area = image_cells_segmented_current_counts[iterator]

                cell_mask = image_cells_segmented_current == current_label    # Mask for elements with current label

                DNA_area_in_cell = image_blue_segmented_current[cell_mask]  # Extract DNA elements in current cell
                DNA_area_in_cell_count = np.count_nonzero(DNA_area_in_cell == 1)  # Count DNA elements with value 1
                
                overlap_area_in_cell = image_overlap_current[cell_mask]
                overlap_area_in_cell_count = np.count_nonzero(overlap_area_in_cell == 1)

                if DNA_area_in_cell_count > 0:
                    overlap_area_percent_of_DNA_area = overlap_area_in_cell_count / DNA_area_in_cell_count * 100
                else:
                    overlap_area_percent_of_DNA_area = np.nan

                DNA_signal_intensity_in_cell_total = np.sum(image_DNA_in_cells_current[cell_mask])

                if current_area > 0:
                    DNA_signal_intensity_in_cell_by_area = DNA_signal_intensity_in_cell_total / current_area
                else:
                    DNA_signal_intensity_in_cell_by_area = np.nan

                cell_df.loc[len(cell_df)] = [int(time_iterator), int(current_label), int(current_area), int(DNA_area_in_cell_count), int(overlap_area_in_cell_count), overlap_area_percent_of_DNA_area, DNA_signal_intensity_in_cell_by_area]
        
        ## for debugging:
        # fig = plt.figure()
        # ax1 = fig.add_subplot(2,1,1)
        # ax1.imshow(image_nuclei_segmented_current, cmap='prism')
        # ax2 = fig.add_subplot(2,1,2)
        # ax2.imshow(image_cells_segmented_current, cmap='prism')
        # #fig.show()
        # plt.show()

        #if verbose mode is active, save intermediate image (labelled cells)
        if verbose:    
            labeled_cells_dir = os.path.join(directory, 'labeleled_cells', ('labeled_cells_' + str(time_iterator) + '.tif'))
            image_cells_segmented_current_out = np.copy(image_cells_segmented_current).astype('uint16')
            tif.imwrite(labeled_cells_dir, image_cells_segmented_current_out)


        print(f'working on timepoint {time_iterator:03} out of {t_max:03} - [x] Cleaning nucleus segmentation image - [x] Segmenting individual cells - [x] Removing {len(edge_total_unique)} cells that touch the edge - [x] Quantifying DNA inside the cells - [x] Characterizing individual cells - [ ] Counting colocalized spots',end='\r')
        
        #assign unique label to each point in the resulting overlap map (gives count)
        image_overlap_labels = np.copy(image_overlap_current)

        image_overlap_labels = skimage.measure.label(image_overlap_labels, background=0)

        if verbose:
            labeled_cells_dir = os.path.join(directory, 'labeled_colocalized_spots', ('labeled_spots_' + str(time_iterator) + '.tif'))
            image_spots_segmented_current_out = np.copy(image_overlap_labels).astype('uint16')
            tif.imwrite(labeled_cells_dir, image_spots_segmented_current_out)

        #use image moments to characterize overlap map (gives area, etc.)
        moments = skimage.measure.regionprops_table(image_overlap_labels, properties=('label','area','centroid'))

        current_label = 0
        max_label = len(moments['label'])

        while current_label < max_label:
            current_area = moments['area'][current_label]
            if image_dimensions == 3:
                current_center_y = int(moments['centroid-0'][current_label])
                current_center_x = int(moments['centroid-1'][current_label])
            elif image_dimensions == 4:
                current_center_z = int(moments['centroid-0'][current_label])
                current_center_y = int(moments['centroid-1'][current_label])
                current_center_x = int(moments['centroid-2'][current_label])
            
            if current_area > overlap_area_cutoff_adjusted:
                (f'Removing overlap label: {current_label:04} area was {current_area}')
            else:     
                if image_dimensions == 3:
                    current_cell_label = int(image_cells_segmented_current[current_center_y,current_center_x])
                elif image_dimensions == 4:
                    current_cell_label = int(image_cells_segmented_current[current_center_z, current_center_y,current_center_x])
                overlap_df.loc[len(overlap_df)] = [time_iterator, current_label, current_cell_label]
            
            current_label += 1

        print(f'working on timepoint {time_iterator:03} out of {t_max:03} - [x] Cleaning nucleus segmentation image - [x] Segmenting individual cells - [x] Removing {len(edge_total_unique)} cells that touch the edge - [x] Quantifying DNA inside the cells - [x] Characterizing individual cells - [x] Counting colocalized spots')

        if count_green_puncta:
            image_cells_segmented_current_mask = np.copy(image_cells_segmented_current)
            image_cells_segmented_current_mask[image_cells_segmented_current_mask > 0] = 1
            image_green_segmented_current_label = np.multiply(image_green_segmented_current, image_cells_segmented_current_mask)

            image_green_segmented_current_label = skimage.measure.label(image_green_segmented_current_label, background = 0)

            moments_green = skimage.measure.regionprops_table(image_green_segmented_current_label, properties=('label','area','centroid'))

            if verbose:    
                labeled_cells_dir = os.path.join(directory, 'labeled_puncta', ('labeled_puncta_' + str(time_iterator) + '.tif'))
                image_green_spots_segmented_current_out = np.copy(image_green_segmented_current_label).astype('uint16')
                tif.imwrite(labeled_cells_dir, image_green_spots_segmented_current_out)

            current_label = 0
            max_label = len(moments_green['label'])

            while current_label < max_label:
                current_area = moments_green['area'][current_label]
                if image_dimensions == 3:
                    current_center_y = int(moments_green['centroid-0'][current_label])
                    current_center_x = int(moments_green['centroid-1'][current_label])
                elif image_dimensions == 4:
                    current_center_z = int(moments_green['centroid-0'][current_label])
                    current_center_y = int(moments_green['centroid-1'][current_label])
                    current_center_x = int(moments_green['centroid-2'][current_label])
                
                if current_area > overlap_area_cutoff_adjusted:
                    (f'Removing green label: {current_label:04} area was {current_area}')
                else:     
                    if image_dimensions == 3:
                        current_cell_label = int(image_cells_segmented_current[current_center_y,current_center_x])
                    elif image_dimensions == 4:
                        current_cell_label = int(image_cells_segmented_current[current_center_z, current_center_y,current_center_x])
                    green_puncta_df.loc[len(green_puncta_df)] = [time_iterator, current_label, current_cell_label]
                
                current_label += 1

        time_iterator += 1

    # Merge dataframes on 'Cell_label' and 'Time' columns
    merged_df = pd.merge(cell_df, overlap_df, on=['Cell_label', 'Time'], how='left')#.fillna(0)
    merged_df = merged_df.dropna(subset=['Overlap_label'])

    # Group by 'Cell_label' and calculate the event count
    event_counts = merged_df.groupby(['Time', 'Cell_label']).size().reset_index(name='Event_count')


    # Add 'Event_count' column to cells_df
    cell_df = pd.merge(cell_df, event_counts, on=['Time','Cell_label'], how='left').fillna(0)

    cell_df['Event_count_per_area'] = cell_df['Event_count'] / cell_df['Cell_area']

    if count_green_puncta:
        merged_df_green = pd.merge(cell_df, green_puncta_df, on=['Cell_label', 'Time'], how='left')
        merged_df_green = merged_df_green.dropna(subset=['Overlap_label'])

        green_puncta_count = merged_df_green.groupby(['Time', 'Cell_label']).size().reset_index(name='Green_puncta_count')
        cell_df = pd.merge(cell_df, green_puncta_count, on=['Time','Cell_label'], how='left').fillna(0)
        cell_df['Green_puncta_count_per_area'] = cell_df['Green_puncta_count'] / cell_df['Cell_area']

    if count_green_puncta:
        result_df = cell_df.groupby('Time').agg({
        'Event_count': ['mean', 'std'],
        'Cell_area': ['mean', 'std', 'count'],
        'Event_count_per_area': ['mean','std'],
        'DNA_area': ['mean','std'],
        'Overlap_area': ['mean','std'],
        'Overlap_area_percent_of_DNA_area': ['mean','std'],
        'DNA_mean_signal_intensity': ['mean','std'],
        'Green_puncta_count': ['mean', 'std'],
        'Green_puncta_count_per_area': ['mean','std']
        })
    else:
        result_df = cell_df.groupby('Time').agg({
        'Event_count': ['mean', 'std'],
        'Cell_area': ['mean', 'std', 'count'],
        'Event_count_per_area': ['mean','std'],
        'DNA_area': ['mean','std'],
        'Overlap_area': ['mean','std'],
        'Overlap_area_percent_of_DNA_area': ['mean','std'],
        'DNA_mean_signal_intensity': ['mean','std']
        })

    if create_histogram == True:
        print('Drawing histograms')
        timepoint_list = np.unique(cell_df['Time'].to_numpy())
        max_cell_count = int(result_df['Cell_area']['count'].max())
        max_event_count = int(cell_df['Event_count'].max())

        for current_timepoint in timepoint_list:
            title_text = f'Event count histogram timepoint {current_timepoint}'
            cell_df_time_subset = cell_df.loc[cell_df['Time'] == current_timepoint]
            cell_df_time_subset_list = cell_df_time_subset['Event_count'].tolist()
            plt.figure().set_figwidth(max_event_count)
            axes = plt.gca()
            axes.set_ylim([0,max_cell_count])
            plt.hist(cell_df_time_subset_list,bins=range(max_event_count),color='blue')
            plt.title(title_text)
            histogram_directory = os.path.join(directory,'histograms', f'hist_{int(current_timepoint):03}.png')
            plt.savefig(histogram_directory)
            plt.close()

    # Load the Excel template file
    xlsx_path = os.path.join(directory, f'{parent_folder_name}_results.xlsx')
    workbook = openpyxl.load_workbook(xlsx_path)

    # Select the right worksheet
    sheet = workbook.worksheets[0]

    # Write DataFrame to Excel file
    for row in result_df.itertuples(index=True):
        sheet.append(list(row))

    sheet = workbook.worksheets[1]
    for row in cell_df.itertuples(index=False):
        sheet.append(list(row))

    sheet = workbook.worksheets[2]
    for row in dna_count_df.itertuples(index=False):
        sheet.append(list(row))

    metadata_list_settings = [['version',version],['nucleus_area_cutoff',nucleus_area_cutoff],['overlap_area_cutoff',overlap_area_cutoff],['z_dimension_scale_factor',z_dimension_scale_factor],['verbose',verbose],['create_histogram',create_histogram],['count_green_puncta',count_green_puncta]]
    sheet = workbook.worksheets[3]
    for row in metadata_list_settings:
        sheet.append(row)

    #save the changes
    workbook.save(xlsx_path)

    return(True)

def main(base_directory):
    # Loop through all subdirectories in the provided directory
    for item in os.listdir(base_directory):
        if os.path.isdir(os.path.join(base_directory,item)):
            sub_directory_path = os.path.join(base_directory,item)
            process_directory(sub_directory_path)

main(data_dir)
