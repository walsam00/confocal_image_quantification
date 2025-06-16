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

import os
import subprocess
import tkinter as tk
from tkinter import filedialog, ttk
import sys
import shutil
import traceback
import configparser

import numpy as np
import pandas as pd
import skimage
import tifffile as tif
import openpyxl
import matplotlib.pyplot as plt


def quantify_tiff_files(cells_seg_tif, nuclei_seg_tif, DNA_tif, DNA_seg_tif, endo_seg_tif, template_dir, verbose, count_green, create_histogram):
    #set default values
    version = '2.00-GUI'
    global nucleus_area_cutoff
    global overlap_area_cutoff
    global z_dimension_scale_factor

    # Copy the template Excel file to the destination path
    base = base_dir.get()
    parent_folder_name = os.path.split(base)[1]
    template_destination_path = os.path.join(base, f'{parent_folder_name}_results.xlsx')
    try:
        shutil.copy(template_dir, template_destination_path)
    except Exception:
        print('failed to copy results template')
        traceback.print_exc()
        return
    
    if verbose == 1:
        try:
            os.mkdir(os.path.join(base, 'labeleled_cells'))
            os.mkdir(os.path.join(base, 'labeled_puncta'))
            os.mkdir(os.path.join(base, 'labeled_colocalized_spots'))
        except Exception:
            print('mkdir labeled_cells failed')
            traceback.print_exc()
            return

    if create_histogram == 1:
        try:
            os.mkdir(os.path.join(base, 'histograms'))
        except Exception:
            print('mkdir histogram failed')
            traceback.print_exc()
            return

    #read input files
    try:
        image_blue_segmented = tif.imread(DNA_seg_tif)
    except Exception:
        print('failed to read segmented DNA image')
        traceback.print_exc()
        return
    try:
        image_blue = tif.imread(DNA_tif)
    except Exception:
        print('failed to read DNA image')
        traceback.print_exc()
        return
    try:
        image_green_segmented = tif.imread(endo_seg_tif)
    except Exception:
        print('failed to read segmented endosomal event image')
        traceback.print_exc()
        return
    try:
        image_nuclei_segmented = tif.imread(nuclei_seg_tif)
    except Exception:
        print('failed to read segmented nuclei image')
        traceback.print_exc()
        return
    try:
        image_cells_segmented = tif.imread(cells_seg_tif)
    except Exception:
        print('failed to read segmented cells image')
        traceback.print_exc()
        return

    #shift the segmentation categories to 0 = background and 1 = foreground
    image_blue_segmented[image_blue_segmented == 1] = 0
    image_blue_segmented[image_blue_segmented == 2] = 1
    image_green_segmented[image_green_segmented == 1] = 0
    image_green_segmented[image_green_segmented == 2] = 1
    image_cells_segmented[image_cells_segmented == 1] = 0
    image_cells_segmented[image_cells_segmented == 2] = 1
    image_nuclei_segmented[image_nuclei_segmented == 1] = 0
    
    #overlap segmentation green & segmentation blue & segmentation cells
    image_overlap = np.multiply(image_green_segmented, image_blue_segmented)
    image_overlap = np.multiply(image_overlap, image_cells_segmented)

    image_DNA_in_cells = np.copy(image_blue)
    image_DNA_in_cells[image_cells_segmented == 0] = 0
    image_DNA_outside_cells = np.copy(image_blue)
    image_DNA_outside_cells[image_cells_segmented != 0] = 0

    #assess image dimensions to get number of timepoints
    image_dimensions = image_overlap.ndim
    if image_dimensions == 3:
        t_max,y_max,x_max = image_overlap.shape
    elif image_dimensions == 4:
        t_max,z_max,y_max,x_max = image_overlap.shape
    time_iterator = 0

    #create dataframes to hold the quantification
    overlap_df = pd.DataFrame(columns=['Time', 'Overlap_label', 'Cell_label'])
    green_puncta_df = pd.DataFrame(columns=['Time', 'Overlap_label', 'Cell_label'])
    cell_df = pd.DataFrame(columns=['Time', 'Cell_label', 'Cell_area', 'DNA_area', 'Overlap_area', 'Overlap_area_percent_of_DNA_area', 'DNA_mean_signal_intensity'])
    dna_count_df = pd.DataFrame(columns=['Time','DNA_in_cells','DNA_outside_cells','area_inside_cells','area_outside_cells','DNA_in_cells_per_area','DNA_outside_cells_per_area'])

    #quantify data timepoint by timepoint
    try:
        while time_iterator < t_max:
            print(f'working on timepoint {time_iterator+1:03} out of {t_max:03} - Cleaning nucleus segmentation image')

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
                
            #label individual nuclei - this assumes the nuclei are never touching, which seems to be a safe assumption
            image_nuclei_segmented_labels = np.copy(image_nuclei_segmented_current)
            image_nuclei_segmented_labels = skimage.measure.label(image_nuclei_segmented_labels, background=0)
            
            #prepare list of nuclei to be removed from the count due to their small size
            nuclei_list_to_clean = []

            #count area of each nucleus, fill in list to be cleaned (exclusion threshold set in settings.ini file)
            moments = skimage.measure.regionprops_table(image_nuclei_segmented_labels, properties=('label','area'))
            
            #check if any nuclei have been found
            if len(moments['label']) == 0:
                raise ValueError(f'Encountered timepoint with no detectable nuclei, cancelling quantification - timepoint {time_iterator}')
            
            current_cell_count = moments['label'][-1]
            for current_label in moments['label']:
                if moments['area'][current_label-1] < nucleus_area_cutoff_adjusted:
                    nuclei_list_to_clean.append(current_label)
                    current_area = moments['area'][current_label-1]
                    current_cell_count -= 1

            #do actual cleaning on the original nucleus segmentation image
            for cleaning_label in nuclei_list_to_clean:
                image_nuclei_segmented_current[image_nuclei_segmented_labels == cleaning_label] = 0

            #label each nucleus in that image
            image_nuclei_segmented_current = skimage.measure.label(image_nuclei_segmented_current, background=0)

            print(f'working on timepoint {time_iterator+1:03} out of {t_max:03} - Segmenting individual cells')

            #perform watershed segmentation using labelled nuclei as troughs, resulting in separated cells even if their borders touch
            image_cells_segmented_current = skimage.segmentation.watershed(image_cells_segmented_current, markers=image_nuclei_segmented_current, mask=(image_cells_segmented_current > 0))

            print(f'working on timepoint {time_iterator+1:03} out of {t_max:03} - Quantifying DNA outside the cells')

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
            elif image_dimensions == 4:
                #edge_1 = image_cells_segmented_current[0,:,:]
                #edge_2 = image_cells_segmented_current[-1,:,:]
                edge_3 = image_cells_segmented_current[:,0,:]
                edge_4 = image_cells_segmented_current[:,-1,:]
                edge_5 = image_cells_segmented_current[:,:,0]
                edge_6 = image_cells_segmented_current[:,:,-1]
                edge_total = np.concatenate((edge_3,edge_4,edge_5,edge_6))#np.concatenate((edge_1,edge_2,edge_3,edge_4,edge_5,edge_6))
            edge_total_unique = np.unique(edge_total)

            print(f'working on timepoint {time_iterator+1:03} out of {t_max:03} - Removing {len(edge_total_unique)} cells that touch the edge')

            for element in edge_total_unique:
                image_cells_segmented_current[image_cells_segmented_current == element] = 0

            print(f'working on timepoint {time_iterator+1:03} out of {t_max:03} - Quantifying DNA inside the cells')

            #count DNA (by summing intesity) inside all cells. Normalize by area
            image_DNA_in_cells_current[image_cells_segmented_current == 0] = 0
            current_DNA_count_in_cells = np.sum(image_DNA_in_cells_current)
            current_area_in_cells = (image_cells_segmented_current > 0).sum()
            current_DNA_count_in_cells_by_area = current_DNA_count_in_cells / current_area_in_cells
            dna_count_df.loc[len(dna_count_df)] = [int(time_iterator), current_DNA_count_in_cells, current_DNA_count_outside_cells, current_area_in_cells, current_area_outside_cells, current_DNA_count_in_cells_by_area, current_DNA_count_outside_cells_by_area]
            
            print(f'working on timepoint {time_iterator+1:03} out of {t_max:03} - Characterizing individual cells')
            
            #get cell count to make cell list
            image_cells_segmented_current_unique, image_cells_segmented_current_counts = np.unique(image_cells_segmented_current, return_counts=True)
            
            #characterize each cell one by one
            max_iterator = len(image_cells_segmented_current_unique)
            if max_iterator == 0:
                raise ValueError(f'Encountered timepoint with no detectable cells, cancelling quantification - timepoint {time_iterator}')
            for iterator in range(max_iterator):
                if iterator > 0:
                    current_label = image_cells_segmented_current_unique[iterator]
                    current_area = image_cells_segmented_current_counts[iterator]

                    #create mask for each cell
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
            
            #if verbose mode is active, save intermediate image (labelled cells)
            if verbose == 1:    
                labeled_cells_dir = os.path.join(base, 'labeleled_cells', ('labeled_cells_' + str(time_iterator) + '.tif'))
                image_cells_segmented_current_out = np.copy(image_cells_segmented_current).astype('uint16')
                tif.imwrite(labeled_cells_dir, image_cells_segmented_current_out)


            print(f'working on timepoint {time_iterator+1:03} out of {t_max:03} Counting colocalized spots',end='\r')
            
            #assign unique label to each point in the overlap map (gives count)
            image_overlap_labels = np.copy(image_overlap_current)
            image_overlap_labels = skimage.measure.label(image_overlap_labels, background=0)

            #if verbose mode is active, save intermediate image (labelled overlapped spots)
            if verbose:
                labeled_cells_dir = os.path.join(base, 'labeled_colocalized_spots', ('labeled_spots_' + str(time_iterator) + '.tif'))
                image_spots_segmented_current_out = np.copy(image_overlap_labels).astype('uint16')
                tif.imwrite(labeled_cells_dir, image_spots_segmented_current_out)

            #use image moments to characterize overlap map (gives area, etc.), if desireable
            moments = skimage.measure.regionprops_table(image_overlap_labels, properties=('label','area','centroid'))

            current_label = 0
            max_label = len(moments['label'])

            #for each colocalized spot, use the image moments derived centroid to get coordinates of each spot's center
            while current_label < max_label:
                current_area = moments['area'][current_label]
                if image_dimensions == 3:
                    current_center_y = int(moments['centroid-0'][current_label])
                    current_center_x = int(moments['centroid-1'][current_label])
                elif image_dimensions == 4:
                    current_center_z = int(moments['centroid-0'][current_label])
                    current_center_y = int(moments['centroid-1'][current_label])
                    current_center_x = int(moments['centroid-2'][current_label])
                
                #Clean the overlap data: remove spots where the overlap is too large (cutoof set in the settings.ini file)
                #once cleaned, assign spots to a cell using center coordinates
                if current_area > overlap_area_cutoff_adjusted:
                    (f'Removing overlap label: {current_label:04} area was {current_area}')
                else:     
                    if image_dimensions == 3:
                        current_cell_label = int(image_cells_segmented_current[current_center_y,current_center_x])
                    elif image_dimensions == 4:
                        current_cell_label = int(image_cells_segmented_current[current_center_z, current_center_y,current_center_x])
                    overlap_df.loc[len(overlap_df)] = [time_iterator, current_label, current_cell_label]
                current_label += 1

            #if the green spots should be counted by themselves also, do so here
            if count_green == 1:
                image_cells_segmented_current_mask = np.copy(image_cells_segmented_current)
                image_cells_segmented_current_mask[image_cells_segmented_current_mask > 0] = 1
                image_green_segmented_current_label = np.multiply(image_green_segmented_current, image_cells_segmented_current_mask)

                image_green_segmented_current_label = skimage.measure.label(image_green_segmented_current_label, background = 0)

                moments_green = skimage.measure.regionprops_table(image_green_segmented_current_label, properties=('label','area','centroid'))

                if verbose:    
                    labeled_cells_dir = os.path.join(base, 'labeled_puncta', ('labeled_puncta_' + str(time_iterator) + '.tif'))
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
        # Normalize by cell area
        cell_df['Event_count_per_area'] = cell_df['Event_count'] / cell_df['Cell_area']

        #add green puncta count if option is enabled
        if count_green == 1:
            merged_df_green = pd.merge(cell_df, green_puncta_df, on=['Cell_label', 'Time'], how='left')
            merged_df_green = merged_df_green.dropna(subset=['Overlap_label'])

            green_puncta_count = merged_df_green.groupby(['Time', 'Cell_label']).size().reset_index(name='Green_puncta_count')
            cell_df = pd.merge(cell_df, green_puncta_count, on=['Time','Cell_label'], how='left').fillna(0)
            cell_df['Green_puncta_count_per_area'] = cell_df['Green_puncta_count'] / cell_df['Cell_area']

        #summarize results by timepoint
        if count_green == 1:
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

        #create overlap count by cell histogram if option is enabled
        if create_histogram == 1:
            print('Drawing histograms')
            timepoint_list = np.unique(cell_df['Time'].to_numpy())
            max_cell_count = int(result_df['Cell_area']['count'].max())
            max_event_count = int(cell_df['Event_count'].max())
            for current_timepoint in timepoint_list:
                title_text = f'Event count histogram timepoint {current_timepoint}'
                cell_df_time_subset = cell_df.loc[cell_df['Time'] == current_timepoint]
                cell_df_time_subset_list = cell_df_time_subset['Event_count'].tolist()
                #max_event_count = int(max(cell_df_time_subset_list))
                plt.figure().set_figwidth(max_event_count)
                axes = plt.gca()
                axes.set_ylim([0,max_cell_count])
                plt.hist(cell_df_time_subset_list,bins=range(max_event_count),color='blue')
                plt.title(title_text)
                histogram_directory = os.path.join(base,'histograms', f'hist_{int(current_timepoint):03}.png')
                plt.savefig(histogram_directory)
                plt.close()

        # Load the Excel file to write results to
        xlsx_path = os.path.join(base, f'{parent_folder_name}_results.xlsx')
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

        metadata_list_settings = [
            ['version',version],
            ['nucleus_area_cutoff',nucleus_area_cutoff],
            ['overlap_area_cutoff',overlap_area_cutoff],
            ['z_dimension_scale_factor',z_dimension_scale_factor],
            ['verbose',verbose],
            ['create_histogram',create_histogram],
            ['count_green_puncta',count_green],
            ['DNA_segmented',DNA_seg_tif],
            ['DNA_raw',DNA_tif],
            ['Endosomal_event_sgmented',endo_seg_tif],
            ['Nuclei_segmented',nuclei_seg_tif],
            ['Cells_segmented',cells_seg_tif]
        ]
        sheet = workbook.worksheets[3]
        for row in metadata_list_settings:
            sheet.append(row)

        #save the changes
        workbook.save(xlsx_path)
    except ValueError as err:
        print('Error in image data: ', err)
        return('red')
    return('green')
        
def run_ilastik_headless(project_name,input_image):
    #run headless ilastik pixel classification
    if project_name.endswith('cells.ilp') or project_name.endswith('cells_3d.ilp'):
        command = [
            ilastik_dir.get(),
            "--headless",
            "--project=" + project_name,
            "--output_format=multipage tiff",
            "--export_source=Simple Segmentation",
            r"--output_filename_format={dataset_dir}\{nickname}_cells_segmented.tif",
            input_image
        ]
    else:    
        command = [
            ilastik_dir.get(),
            "--headless",
            "--project=" + project_name,
            "--output_format=multipage tiff",
            "--export_source=Simple Segmentation",
            r"--output_filename_format={dataset_dir}\{nickname}_segmented.tif",
            input_image
        ]

    try:
        subprocess.run(command, check=True)
        print("Ilastik processing completed successfully!")
        print(f'command: {command}')
    except subprocess.CalledProcessError as e:
        print(f"Error running Ilastik: {e}")

def select_directory():
    directory = filedialog.askdirectory()
    if directory:
        base_dir.set(directory)
        update_initial_files()

def update_initial_files():
    #pre-populate all the filepaths for which actual files are found in the specified directory
    ilastik_file_list_2D = ["pixel_classificaiton_green_cells.ilp", "pixel_classificaiton_nuclei.ilp", "pixel_classificaiton_blue.ilp", "pixel_classificaiton_green.ilp"]
    ilastik_file_list_3D = ["pixel_classificaiton_green_cells_3d.ilp", "pixel_classificaiton_nuclei_3d.ilp", "pixel_classificaiton_blue_3d.ilp", "pixel_classificaiton_green_3d.ilp"]
    base = base_dir.get()
    for i in range(len(output_files)):
        output_files[i].set('')
    if base:
        output_template_file.set(os.path.join(os.path.split(__file__)[0],'template_per_cell_data.xlsx'))
        for i, filename in enumerate(["green.tif", "nuclei.tif", "blue.tif", "green.tif"]):
            matching_files = [f for f in os.listdir(base) if f.endswith(filename)]
            if matching_files:
                file_path = os.path.abspath(os.path.join(base, matching_files[0]))  # Use the first match
                file_path_ilastik = os.path.join(os.path.split(__file__)[0], ilastik_file_list_2D[i])
                input_files[i].set(file_path)
                input_ilastik_files[i].set(file_path_ilastik)
                process_flags[i].set(1)  # Enable processing for existing files by default
                if i==0:
                    if os.path.isfile(input_files[i].get().replace(".tif", "_cells_segmented.tiff")):
                        output_files[i].set(input_files[i].get().replace(".tif", "_cells_segmented.tiff"))
                else:
                    if os.path.isfile(input_files[i].get().replace(".tif", "_segmented.tiff")):
                        output_files[i].set(input_files[i].get().replace(".tif", "_segmented.tiff"))
            else:
                filenam3 = filename.replace('e','3')
                matching_files = [f for f in os.listdir(base) if f.endswith(filenam3)]
                if matching_files:
                    file_path = os.path.join(base, matching_files[0])  # Use the first match
                    file_path_ilastik = os.path.join(os.path.split(__file__)[0], ilastik_file_list_3D[i])      #os.path.abspath(os.path.join(base, "..", ilastik_file_list_3D[i]))
                    input_files[i].set(file_path)
                    input_ilastik_files[i].set(file_path_ilastik)
                    process_flags[i].set(1)  # Enable processing for existing files
                    if i==0:
                        if os.path.isfile(input_files[i].get().replace(".tif", "_cells_segmented.tiff")):
                            output_files[i].set(input_files[i].get().replace(".tif", "_cells_segmented.tiff"))
                    else:
                        if os.path.isfile(input_files[i].get().replace(".tif", "_segmented.tiff")):
                            output_files[i].set(input_files[i].get().replace(".tif", "_segmented.tiff"))
                else:
                    input_files[i].set("")
                    process_flags[i].set(0)  # Disable processing if file does not exist
    segmentation_button.config(bg='white')
    quantification_button.config(bg='white')
        

def select_file(var):
    filename = filedialog.askopenfilename()
    if filename:
        var.set(filename)

def select_ilastik_executable(var):
    filename = filedialog.askopenfilename()
    if filename:
        ilastik_dir.set(filename)
        config = configparser.ConfigParser()
        config_path = 'GUI_settings.ini'  #os.path.join( os.path.split(base_dir.get())[1],'GUI_settings.ini')
        config.read(config_path)
        config.set('Settings', 'ilastik_executable', filename)

        with open(config_path, 'w') as configfile:
            config.write(configfile)

def process_step_one():
    #initialize image segmentation for all enabled files
    for i in range(4):
        if process_flags[i].get() and input_files[i].get():
            print(f"Processing: {input_files[i].get()}")
            run_ilastik_headless(input_ilastik_files[i].get(), input_files[i].get())

            if i==0:
                if os.path.isfile(input_files[i].get().replace(".tif", "_cells_segmented.tiff")):
                    output_files[i].set(input_files[i].get().replace(".tif", "_cells_segmented.tiff"))
            else:
                if os.path.isfile(input_files[i].get().replace(".tif", "_segmented.tiff")):
                    output_files[i].set(input_files[i].get().replace(".tif", "_segmented.tiff"))
    segmentation_button.config(bg='green')
                

def process_step_two():
    #initialize image quantification
    quantification_button.config(bg='white')
    if os.path.isfile(output_files[0].get()):
        if os.path.isfile(output_files[1].get()):
            if os.path.isfile(output_files[2].get()):
                if os.path.isfile(output_files[3].get()):
                    if os.path.isfile(output_template_file.get()):
                        outcome = quantify_tiff_files(output_files[0].get(), output_files[1].get(), input_files[2].get(), output_files[2].get(), output_files[3].get(), output_template_file.get(), quant_flags[0].get(), quant_flags[1].get(), quant_flags[2].get())
                    else:
                        print('Output Excel template not found')
                else:
                    print('Endosomal event segmentation file not found')
            else:
                print('DNA segmentation file not found')
        else:
            print('Nuclei segmentation file not found')
    else:
        print('Cell segmentation file not found')
    quantification_button.config(bg=outcome)


def show_tooltip(event, text):
    #create tooltip on mouse hover to display long filenames
    tooltip = tk.Toplevel(root)
    tooltip.wm_overrideredirect(True)
    tooltip.wm_geometry(f"+{event.x_root+10}+{event.y_root+10}")
    label = tk.Label(tooltip, text=text, background="yellow", relief="solid", borderwidth=1)
    label.pack()
    event.widget.tooltip = tooltip

def hide_tooltip(event):
    #hide tooltip when mouse leaves relevant widget
    if hasattr(event.widget, "tooltip"):
        event.widget.tooltip.destroy()
        event.widget.tooltip = None

#------------------------------------------------------------------------------------------------#
# GUI creation section, organized by frame
#------------------------------------------------------------------------------------------------#
root = tk.Tk()
root.title("Confocal image quantification")
icon = tk.PhotoImage(file='icon.png')
root.iconphoto(False, icon)

base_dir = tk.StringVar()

root.columnconfigure(0, minsize=800, weight=1)

#------------------------------------------------------------------------------------------------#

frm_input_selection = tk.Frame(master=root, padx=10)
frm_input_selection.columnconfigure(1, weight=1)

tk.Button(frm_input_selection, text="Select Base Directory", command=select_directory).grid(row=0, column=0, columnspan=1, sticky='ew')
tk.Entry(frm_input_selection, textvariable=base_dir, width=50).grid(row=0, column=1, columnspan=3, sticky='enw')
ilastik_dir = tk.StringVar()
tk.Button(frm_input_selection, text="Ilastik executable", command=lambda var=ilastik_dir: select_ilastik_executable(var)).grid(row=1, column=0, columnspan=1,sticky='ew')
tk.Entry(frm_input_selection, textvariable=ilastik_dir, width=50).grid(row=1, column=1, columnspan=3, sticky='enw')

#------------------------------------------------------------------------------------------------#

input_files = [tk.StringVar() for _ in range(4)]
input_ilastik_files = [tk.StringVar() for _ in range(4)]
output_files = [tk.StringVar() for _ in range(4)]
process_flags = [tk.IntVar(value=1) for _ in range(4)]  # Default all to 1 (checked)
quant_flags = [tk.IntVar(value=0) for _ in range(3)]
output_template_file = tk.StringVar()
labels = ["Cells channel", "Nuclei channel", "DNA channel", "Endosomal event channel"]

### Input Parameters ###
config = configparser.ConfigParser()
config_path = 'GUI_settings.ini'
config.read(config_path)
nucleus_area_cutoff = float(config['Settings']['Nucleus_cutoff_sidelength'])            #sidelength of the square (for 2D) or the cube (in 3D) that defines the minimal acceptable area (in pixels) of a nucleus. Nuclei smaller than this are removed before analysis
overlap_area_cutoff = float(config['Settings']['Overlap_cutoff_sidelength'])            #sidelength of the square (for 2D) or the cube (in 3D) that defines the maximum acceptable sized of colocalized areas (in pixels) of blue and green. Spots larger than this are removed before analysis
z_dimension_scale_factor = float(config['Settings']['z_dimension_scale_factor'])        #how much smaller the scale of the x/y-axis is (for z-stacks) than of the z-axis -> is used to adequately calculate the overlap/nucleus cutoff values for 3D
ilastik_dir.set(config['Settings']['ilastik_executable'])

frm_segmentation = tk.Frame(master=root, padx=10)
frm_segmentation.columnconfigure(2, weight=1)
frm_segmentation.columnconfigure(1, minsize=200)

for i, label in enumerate(labels):
    tk.Checkbutton(frm_segmentation, variable=process_flags[i]).grid(row=2*i, column=0, sticky='w')
    tk.Label(frm_segmentation, text=label+' tif', anchor='w').grid(row=2*i, column=1, sticky="we")
    entry = tk.Entry(frm_segmentation, textvariable=input_files[i])
    entry.grid(row=2*i, column=2, sticky="ew")
    entry.bind("<Enter>", lambda e, var=input_files[i]: show_tooltip(e, var.get()))
    entry.bind("<Leave>", hide_tooltip)
    tk.Button(frm_segmentation, text="Select", command=lambda var=input_files[i]: select_file(var)).grid(row=2*i, column=3, sticky='enw')
    tk.Label(frm_segmentation, text=label + " ilp").grid(row=2*i+1, column=1, sticky="w")
    entry2 = tk.Entry(frm_segmentation, textvariable=input_ilastik_files[i])
    entry2.grid(row=2*i+1, column=2, sticky="we")
    entry2.bind("<Enter>", lambda e, var=input_ilastik_files[i]: show_tooltip(e, var.get()))
    entry2.bind("<Leave>", hide_tooltip)
    tk.Button(frm_segmentation, text="Select", command=lambda var=input_ilastik_files[i]: select_file(var)).grid(row=2*i+1, column=3, sticky='enw')

segmentation_button = tk.Button(frm_segmentation, text="Run segmentation", command=process_step_one, font='TkDefaultFont 10 bold', pady=5)
segmentation_button.grid(row=2*len(labels)+2, column=0, columnspan=4, sticky='w')

#------------------------------------------------------------------------------------------------#

frm_quantification = tk.Frame(master=root, padx=10)
frm_quantification.columnconfigure(2, weight=1)
frm_quantification.columnconfigure(0, minsize=30)
frm_quantification.columnconfigure(1, minsize=200)

for i, label in enumerate(labels):
    tk.Label(frm_quantification, text=label + " Processed", anchor='w').grid(row=i, column=1, sticky="ew")
    entry = tk.Entry(frm_quantification, textvariable=output_files[i])
    entry.grid(row=i, column=2, sticky="ew")
    entry.bind("<Enter>", lambda e, var=output_files[i]: show_tooltip(e, var.get()))
    entry.bind("<Leave>", hide_tooltip)
    tk.Button(frm_quantification, text="Select", command=lambda var=output_files[i]: select_file(var)).grid(row=i, column=3, sticky='enw')

tk.Label(frm_quantification, text='Result template', anchor='w').grid(row=len(labels)+1, column=1, sticky="w")
entry = tk.Entry(frm_quantification, textvariable=output_template_file)
entry.grid(row=len(labels)+1, column=2, sticky="ew")
entry.bind("<Enter>", lambda e, var=output_template_file: show_tooltip(e, var.get()))
entry.bind("<Leave>", hide_tooltip)
tk.Button(frm_quantification, text="Select", command=lambda var=output_template_file: select_file(var)).grid(row=len(labels)+1, column=3, sticky='enw')

tk.Label(frm_quantification, text='Save intermediates', anchor='w').grid(row=len(labels)+2, column=1, sticky='w')
tk.Checkbutton(frm_quantification, variable=quant_flags[0]).grid(row=len(labels)+2, column=2, sticky='w')
tk.Label(frm_quantification, text='Count green puncta individually', justify='left').grid(row=len(labels)+3, column=1, sticky='w')
tk.Checkbutton(frm_quantification, variable=quant_flags[1],offvalue=0, onvalue=1).grid(row=len(labels)+3, column=2, sticky='w')
tk.Label(frm_quantification, text='Create histograms', justify='left').grid(row=len(labels)+4, column=1, sticky='w')
tk.Checkbutton(frm_quantification, variable=quant_flags[2],offvalue=0, onvalue=1).grid(row=len(labels)+4, column=2, sticky='w')

quantification_button = tk.Button(frm_quantification, text="Run quantification", command=process_step_two, font='TkDefaultFont 10 bold', pady=5)
quantification_button.grid(row=len(labels)+5, column=0, columnspan=4, sticky='w')

#------------------------------------------------------------------------------------------------#

frm_input_selection.grid(row=0, column=0, sticky='nwe')

frm_segmentation.grid(row=1, column=0, sticky='nwe')

frm_quantification.grid(row=2, column=0, sticky='nwe')

#------------------------------------------------------------------------------------------------#

root.mainloop()
