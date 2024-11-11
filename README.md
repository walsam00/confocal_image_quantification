# confocal_image_quantification

This software takes time-resolved confocal microscopy images, separated by channel and provided as *.tif files, and detects the amount and characteristics of spots colocalized in two channels. 
The software does this quantification on a per-cell basis and uses cell nuclei and a cell mask it creates from one of the two channel images containing the spots of interest.
Time series are expected as input data. [Ilastik](https://www.ilastik.org/) is used to segment the confocal images based on user-provided training. The software handles 2D as well as 3D (z-stack) confocal images and is intended to run on Windows 10/11.

This software was developed for a project quantifying endosomal processes following DNA-LNP transfection of specific cell lines. Read the paper [here](insert DOI link once published). The software can be used as is using the pre-trained Ilastik image segmentation available [here](link to Ilastik project files once published) to replicate the published data or to apply to more image data of the same type. The software can be adapted without altering the code to any image quantification task that uses the same features (quantifying the overlap of puncta contained in of two signal channels of interest, one additionally containing a weaker signal for the cells themselves, as well as a cell nucleus image channel) by training new Ilastik project files for each image channel.

## How to use:

1. Save confocal_image_quantification folder contained in this repository including its contents on your machine.
2. Download pre-trained Ilastik segmentation project files for each image channel [here](link to Ilastik project files once published) or create your own. Save these project files with the correct file-names in the same location as the .bat file (see further on).
3. From your image data to be analyzed, manually create separate .tif files out of the channels of interest. The software was developed for microscopy images from an Olympus confocal microscope, which produces .oir files. [FIJI](https://imagej.net/software/fiji/downloads) is a good program to split out the channels of interest. The [Olympus FIJI plugin](https://imagej.net/formats/olympus) is the best way to import .oir files into FIJI. To use the software as intended, the following image channels are needed as separate .tif-files: 
		a) Labeled DNA signal ('blue' signal)
		b) Signal for events of interest (endosomal escape or recycling) which includes a weaker general staining of the whole cell ('green' signal)
		c) Labeled nuclei signal ('nuclei' signal)
4. Create a new subfolder for those files in the 'automated_colocalization_quantification' Folder.
5. Put all three files into that folder. They need to contain keywords in their filname: 'blue' 'green' 'nuclei'. Z-stacks need to be identified by using the keywords 'blu3', 'gr33n' and 'nucl3i' instead.
6. Run the 'Analysis_per_cell.bat' script. 
		The script automatically applies the 'pixel_classification_blue.ilp' and the 'pixel_classification_green.ilp', 'pixel_classification_nuclei.ilp' and 'pixel_classification_green_cells.ilp' Ilastik project files to segment all unsegmented tif files.
		The batch file then runs the Python script 'Colocalization_quantification_per_cell.py' which overlaps the three masks (blue dots, green dots, cells), creates a 'results_per_cell.xlsx' from the 'template_per_cell_data.xlsx' Excel sheet. The script finally quantifies the result per cell and reports the data by copying it into an Excel template file.

If desired, the following settings can be changed in the settings.ini file (open in Notepad):

+ nucleus_area_cutoff : sidelength of the square (for 2D) or the cube (in 3D) that defines the minimal acceptable area (in pixels) of a nucleus. Nuclei smaller than this are removed before analysis
+ overlap_area_cutoff : sidelength of the square (for 2D) or the cube (in 3D) that defines the maximum acceptable sized of colocalized areas (in pixels) of blue and green. Spots larger than this are removed before analysis
+ z_dimension_scale_factor : factor by how much smaller the scale of the x/y-axis is (for z-stacks) compared to the z-axis -> is used to adequately calculate the overlap/nucleus cutoff values for 3D
+ verbose = if set to 'True', creates subdirectory and saves overlap .tif, green_puncta .tif and labelled cell .tif images. Useful for validating results/debugging/visualization
+ create_histogram : if set to 'True', creates subdirectory and saves histograms of events per cell data per timepoint as png images. Useful for validation of the quantification as well as general science
+ count_green_puncta : if set to 'True', counts green signal puncta in addition to overlap (green AND blue) puncta.

### Prerequisites:
+ Ilastik installation [tested with version 1.4.0]
+ Python 3.x installation
+ Python packages listed further on
+ Ilastik pixel classification project files, trained on representative sample of the data. One for blue spots, one for green spots, one for whole cells, one for nuclei.
+ Correct assignment of file path to the Ilastik executable in the 'Analysis_per_cell.bat' file (use notepad to edit the file).
	
### Important:

+ The script works with relative file paths. The location of the ilastik project files, the '.bat' and '.py' files as well as the subfolders relative to each other cannot be changed, otherwise the software breaks.
+ The location of the ilastik executable that will do the actual segmentation needs to be correctly specified in the 'Analysis_per_cell.bat' file (use notepad to edit the file).
+ The script scans all subdirectories and only segments/quantifies ones it has not already segmented/quantified. In theory you can leave all quantified folders as they are and just keep creating new ones.

### Dependencies:

The Python script needs the following libraries:
+ numpy (tested with version 1.24.3)
+ pandas (tested with version 2.2.1)
+ scikit-image (tested with version 0.22.0)
+ tifffile (tested with version 2024.2.12)
+ openpyxl (tested with version 3.1.2)
+ matplotlib (tested with version 3.8.3)
+ (os, shutil, configparser and glob are also used, but those are contained in the base python install)

Download python 3.* itself via Windows app store (tested with version 3.11.9)

install packages using Windows command prompt (consider creating a dedicated virtual environment):
+ python -m pip install --upgrade pip
+ python -m pip install numpy
+ python -m pip install pandas
+ python -m pip install scikit-image
+ python -m pip install tifffile
+ python -m pip install openpyxl
+ python -m pip install matplotlib

## Required files and folder structure:

Files needed (file- and directory names need to match exactly, unless specified):
|Folder/File|Description|
|---|---|
|**automated_colocalization_quantification** | Parent folder, can be placed anywhere, can be renamed|
|Colocalization_quantification_per_cell.py | Quantification python script|
|Analysis_per_cell.bat								| Execute this to run the software. Windows batch script, first runs Ilastik segmentation, then initiates above python script.|
|settings.ini									| File containing settings for the python script that quantifies the segmented image data. The user can change these settings by editing this file in a word processor (Notepad for example). The format needs to be kept exactly as it is, only the individual values for the settings can be adapted (e.g. changing 'False' to 'True' or changing the numerical value of a parameter).|
|pixel_classification_blue.ilp							| Ilastik pixel classification project file for blue spots, needs to be trained on representative sample of all data to be segmented. Do this once and don't change anymore to ensure consistent segmentation across all files. The filename of this ilastik project file needs to be exactly as written here.|
|pixel_classification_green.ilp							| Ilastik pixel classification project file for green spots, needs to be trained on representative sample of all data to be segmented. Do this once and don't change anymore to ensure consistent segmentation across all files. The filename of this ilastik project file needs to be exactly as written here.|
|pixel_classification_green_cells.ilp						| Ilastik pixel classification project file for whole cells, needs to be trained on representative sample of all data to be segmented. Do this once and don't change anymore to ensure consistent segmentation across all files. The filename of this ilastik project file needs to be exactly as written here.|
|pixel_classification_nuclei.ilp							| Ilastik pixel classification project file for nuclei, needs to be trained on representative sample of all data to be segmented. Do this once and don't change anymore to ensure consistent segmentation across all files. The filename of this ilastik project file needs to be exactly as written here.|
|pixel_classification_blue_3d.ilp						| Ilastik pixel classification project file for blue spots in z-stack images, needs to be trained on representative sample of all z-stack data to be segmented. Do this once and don't change anymore to ensure consistent segmentation across all files. The filename of this ilastik project file needs to be exactly as written here.|
|pixel_classification_green_3d.ilp						| Ilastik pixel classification project file for green spots in z-stack images, needs to be trained on representative sample of all z-stack data to be segmented. Do this once and don't change anymore to ensure consistent segmentation across all files. The filename of this ilastik project file needs to be exactly as written here.|
|pixel_classification_green_cells_3d.ilp						| Ilastik pixel classification project file for whole cells in z-stack images, needs to be trained on representative sample of all z-stack data to be segmented. Do this once and don't change anymore to ensure consistent segmentation across all files. The filename of this ilastik project file needs to be exactly as written here.|
|pixel_classification_nuclei_3d.ilp						| Ilastik pixel classification project file for nuclei in z-stack images, needs to be trained on representative sample of all z-stack data to be segmented. Do this once and don't change anymore to ensure consistent segmentation across all files. The filename of this ilastik project file needs to be exactly as written here.|
|template_per_cell_data.xlsx							| Excel template that is copied and filled in for each set of data that is quantified. Can be adapted and chaged (e.g. the graphs, column titles, etc.), but the style and placement of the data that is added will always be the same.|
|**your data folder pre analysis**| Add this folder yourself: It should contain data to be analyzed. The script can handle any number of these folders at once, just add one folder for each time-series. Folder name can be changed. Must contain the files specified in the following lines, and nothing else:|
|filename_blue.tif							| Blue channel image data as a tif file. Axis order needs to be [t,y,x] or [t,z,y,x] (images must be a time-series). The base filename is ideally the same for all three image files. Combine the base filename and the appropriate keyword ('blue'/'green'/'nuclei' OR 'blu3'/'gr33n'/'nucl3i'). The base filename can be anything, but cannot contain the key words (i.e. the base filename shouldn't contain 'blue'/'green'/'nuclei'/'blu3'/'gr33n'/'nucl3i').|
|filename_green.tif							| Green channel image data as a tif file. Axis order needs to be [t,y,x] or [t,z,y,x] (images must be a time-series). The base filename is ideally the same for all three image files. Combine the base filename and the appropriate keyword ('blue'/'green'/'nuclei' OR 'blu3'/'gr33n'/'nucl3i'). The base filename can be anything, but cannot contain the key words (i.e. the base filename shouldn't contain 'blue'/'green'/'nuclei'/'blu3'/'gr33n'/'nucl3i').|
|filename_nuclei.tif							| Nucleus stain channel image data as a tif file. Axis order needs to be [t,y,x] or [t,z,y,x] (images must be a time-series). The base filename is ideally the same for all three image files. Combine the base filename and the appropriate keyword ('blue'/'green'/'nuclei' OR 'blu3'/'gr33n'/'nucl3i'). The base filename can be anything, but cannot contain the key words (i.e. the base filename shouldn't contain 'blue'/'green'/'nuclei'/'blu3'/'gr33n'/'nucl3i').|
|**your data folder post analysis**						| This is what the data folders look like after the script successfully ran. If the script is ran again with new subfolders in the directory, fully analyzed folders like this one are ignored, and thus can be left where they are, if desired.|
|filename_blue.tif							| Blue channel image data as a tif file. Unchanged post-analysis.|
|filename_blue_segmented.tiff						| Blue channel image data in segmented form (background = 1, spots of interest = 2) as a tiff file.|
|filename_green.tif							| Green channel image data as a tif file. Unchanged post-analysis.|
|filename_green_cells_segmented.tiff					| Green channel image data in segmented form (background = 1, cells = 2) as a tiff file.|
|filename_green_segmented.tiff						| Green channel image data in segmented form (background = 1, spots of interest = 2) as a tiff file.|
|filename_nuclei.tif							| Nucleus stain channel image data as a tif file. Unchanged post-analysis.|
|filename_nuclei_segmented.tiff						| Green channel image data in segmented form (background = 1, nuclei = 2) as a tiff file.|
|results_per_cell.xlsx							| Excel file containing a summary of the analysis as well as detailed data on the second worksheet. Filled copy of template_per_cell_data.xlsx.|
|**your data folder post analysis cleaned**				| If you want to minimize the storage footprint of your data, you can clean up your finished folders like this (script will recognize it as fully analyzed in this state):|
|filename_blue_segmented.tiff						| Optionally keep this segmented data file. Useful if other analyses are planned with this data.|
|filename_green_cells_segmented.tiff					| Optionally keep this segmented data file. Useful if other analyses are planned with this data.|
|filename_green_segmented.tiff						| Optionally keep this segmented data file. Useful if other analyses are planned with this data.|
|filename_nuclei_segmented.tiff					| Optionally keep this segmented data file. Useful if other analyses are planned with this data.|
|results_per_cell.xlsx      | Keep the results Excel file.|
