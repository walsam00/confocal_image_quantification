# confocal_image_quantification
This software takes time-resolved confocal microscopy images, separated by channel and provided as *.tif files, and detects the amount and characteristics of spots co-localized in two of the image channels. 
The software does this quantification on a per-cell basis and uses nuclei and a cell mask to detect each cell. The cell nucleus mask is created from the nucleus confocal image channel; the cell mask is created from a combined cell/endosomal event puncta image channel.
Time series are expected as input data. The [Ilastik](https://www.ilastik.org/) software is used to segment the confocal images based on user-provided training. The software handles 2D as well as 3D (z-stack) confocal images and is intended to run on Windows 10/11. Each segmentation mask requires an Ilastik pixel classification project file, which has been trained by hand on a representative subset of the image data to which the software will be applied. 3D image datasets require separate Ilastik project files from 2D image datasets.

This software was developed for a project quantifying endosomal processes following DNA-LNP transfection of specific cell lines. Read the paper [here](insert DOI link once published). The software can be used as is using the pre-trained Ilastik image segmentation available [here](link to Ilastik project files once published) to replicate the published data or to apply to more image data of the same type. The software can be adapted without altering the code to any image quantification task that requires the same workflow on the same image features (i.e. quantifying the overlap of puncta contained in of two signal channels of interest, with one channel additionally containing a weaker signal for the cells themselves, as well as a separate cell nucleus image channel). The analysis workflow of the software can be applied to different data by training new Ilastik project files for each image channel.

![](/assets/labeled_cells.png)
![](/assets/Event_graph.png)
![](/assets/labeled_cells_3D.png)
![](/assets/event_graph_3D.png)

# How to install:
1. Download the confocal_image_quantification folder contained in this repository including its contents on your machine.
2. Download [Ilastik](https://www.ilastik.org/)
3. Download pre-trained Ilastik segmentation project files for each image channel [here](link to Ilastik project files once published) or create your own. By default, the software expects to find these project files with the correct file-names in the same location as the .bat file used to run the software (details - see further on).
4. Download Python and the required dependencies (see ‘Dependencies’ section)

# Image Data Preparation and Software Usage

## 1. Prepare Your Image Data
- Separate the channels of interest into individual `.tif` files. Use [FIJI](https://imagej.net/software/fiji/downloads) with the [Olympus plugin](https://imagej.net/formats/olympus) to split the channels from Olympus `.oir` files.
- The required channels are:
  - **DNA signal ('blue')**
  - **Events of interest (e.g., endosomal escape) with general cell staining ('green')**
  - **Nuclei signal ('nuclei')**
- Save these `.tif` files in the same folder with the following names:
  - **'blue.tif'**, **'green.tif'**, **'nuclei.tif'** (for single images)
  - **'blu3.tif'**, **'gr33n.tif'**, **'nucl3i.tif'** (for Z-stacks)

## 2. Prepare Ilastik Project Files
- The Ilastik project files should also be in the same folder as the software. Use these names:
  - **'pixel_classification_blue.ilp'**, **'pixel_classification_green.ilp'**, **'pixel_classification_green_cells.tif'**, **'pixel_classification_nuclei.ilp'**
- For Z-stacks, use:
  - **'pixel_classification_blue_3D.ilp'**, **'pixel_classification_green_3D.ilp'**, **'pixel_classification_green_cells_3D.tif'**, **'pixel_classification_nuclei_3D.ilp'**

## 3. Run the Software
- Execute **‘run_confocal_image_quantification_per_cell.bat’**.
- The first time, specify the path to **Ilastik.exe** using the **‘Ilastik executable’** button.

## 4. Select the Base Directory
- Choose the folder where your image files are stored.
- The software will automatically detect the required files if they are correctly named. If not, manually select them.

## 5. Run Image Segmentation
- Click the **‘Run segmentation’** button to process the images with Ilastik.
- After segmentation, the software automatically populates the quantification section.

## 6. Quantification Setup
- The Excel template for results is automatically set. Change it with the **‘Select’** button if needed.
- Enable **‘Save intermediates’**, **‘Count green puncta individually’**, or **‘Create histograms’** for more details or troubleshooting.

## 7. Run Quantification
- Click **‘Run quantification’** to analyze the images. The software needs the segmented image channels and raw DNA channel for analysis.
- Results will be saved in the selected directory.

If desired, the following software settings can be changed in the GUI_settings.ini file (open in Notepad):

+ nucleus_area_cutoff : sidelength of the square (for 2D) or the cube (in 3D) that defines the minimal acceptable area (in pixels) of a nucleus. Nuclei smaller than this are removed before analysis
+ overlap_area_cutoff : sidelength of the square (for 2D) or the cube (in 3D) that defines the maximum acceptable sized of colocalized areas (in pixels) of blue and green. Spots larger than this are removed before analysis
+ z_dimension_scale_factor : factor by how much smaller the scale of the x/y-axis is (for z-stacks) compared to the z-axis -> is used to adequately calculate the overlap/nucleus cutoff values for 3D
+ Ilastik_executable : location of the ‘ilastik.exe’ file. This location can be specified in the software directly.

# Dependencies:
+ Ilastik installation [tested with version 1.4.0]
+ Ilastik pixel classification project files, trained on representative sample of the data. One for blue puncta, one for green puncta, one for whole cells, one for nuclei.
+ Python 3.x installation [tested with version 3.11.9]
The Python script needs the following libraries:
+ numpy (tested with version 1.24.3)
+ pandas (tested with version 2.2.1)
+ scikit-image (tested with version 0.22.0)
+ tifffile (tested with version 2024.2.12)
+ openpyxl (tested with version 3.1.2)
+ matplotlib (tested with version 3.8.3)
+ (os, shutil, configparser, subprocess, sys, traceback, and tkinter are also used, but those are contained in the base python install)

Download python 3.* itself via Windows app store (tested with version 3.11.9)

Install packages using Windows command prompt (consider creating a dedicated virtual environment):
+ python -m pip install --upgrade pip
+ python -m pip install numpy
+ python -m pip install pandas
+ python -m pip install scikit-image
+ python -m pip install tifffile
+ python -m pip install openpyxl
+ python -m pip install matplotlib
