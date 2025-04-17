# Software Overview

This software analyzes time-resolved confocal microscopy images that are separated by channel and provided as `.tif` files. It detects and quantifies the amount and characteristics of spots that are co-localized in two specific image channels. The quantification is performed on a **per-cell** basis, using nuclei and cell masks to identify each cell in the image. Read the related [research article](https://doi.org/10.1016/j.jconrel.2025.113709) for detailed information.

- **Green channel**: Contains cells and endosomal event as puncta
- **Blue channel**: Contains DNA signal that colocalizes with endosomal event puncta
- **Nucleus channel**: Contains cell nuclei, used for separation of individual cells

### Input Data:
- The software expects **time series** of images as input.
- The **Ilastik** software is used for image segmentation based on user-provided training data.
- It works with both **2D** and **3D (z-stack)** confocal images.

### Segmentation Process:
- Each segmentation mask requires an **Ilastik pixel classification project file**. This file is generated through training on a representative subset of the image data to which the software will be applied.
- For **3D** image datasets, separate Ilastik project files are needed, different from the files used for 2D datasets.

### System Requirements:
- The software is designed to run on **Windows 10/11**.

# What the software does

The software quantifies the overlap of the signal in the 'blue' channel with the puncta in the 'green' image channel on a per-cell basis over a time-series of images:
- The three image channels are segmented into four masks by user-trained image segmentation in Ilastik
- The **cell area mask** and the **cell nucleus mask** are combined to separate individual cells using watershed segmentation
  - Nuclei that are below a set size threshold are ignored (adjustable in GUI_settings.ini)
  - Cells that are cut off, i.e. touch the edge of the image (2D - any edge; 3D - any of the sides of the image volume, not top or bottom) are excluded from the analysis
- The **'blue' signal mask** and the **'green puncta' signal mask** are overlapped
- This overlap is quantified on a per-cell basis for each timepoint
    - Individual overlapping puncta are excluded from analysis if they are above a set area/volume cutoff (adjustable in GUI_settings.ini)
- The results are summarized and graphed in an Excel-report and include:
    - Mean number of overlap events per cell over time
    - Mean overlap events per cell and cell area(2D)/volume(3D) over time
    - Mean 'blue' signal intensity per cell over time
    - Mean overlap area/volume per cell over time
    - Mean overlap area/volume in relation to 'blue' signal area/volume (%) per cell over time
    - Mean number of 'green' signal puncta per cell over time
    - Mean number of 'green' signal puncta per cell and cell area over time
    - Mean cell area/volume over time
    - Mean overall intracellular 'blue' signal intensity over time
    - Mean overall extracellular 'blue' signal intensity over time
- The results additionally contain the raw data for each cell at each timepoint

# Use Case and Adaptability

This software was originally developed for a project quantifying endosomal processes following **DNA-LNP transfection** of specific cell lines. The full paper can be found [here](https://doi.org/10.1016/j.jconrel.2025.113709).

You can use the software as-is with the **pre-trained Ilastik segmentation files** available [here](https://doi.org/10.5281/zenodo.14537397) (you need the *.ilp files contained in that repository) to replicate the published results or apply it to similar types of image data. 

While the software is tailored for this specific use case, it can be adapted for other image quantification tasks. The workflow can be used for any task that involves quantifying the overlap of puncta in two signal channels, where one channel also contains a weaker signal representing the cells, along with a separate channel for the cell nucleus. 

To apply the software to different datasets, you can train new Ilastik project files for each image channel.

# Example 2D segmented cell mask (single timepoint of a time-series): #

![](/assets/labeled_cells.png)

# Corresponding overlap event quantification over time: #

![](/assets/Event_graph.png)

Note the size of the error bars. The biological response varies greatly cell-by-cell, which was discovered thanks to the large amount of cells quantified by this software.

# Example 3D segmented cell mask (single timepoint of a time-series): #

![](/assets/labeled_cells_3D.png)

# Corresponding overlap event quantification over time: #

![](/assets/event_graph_3D.png)

Again, note the size of the error bars showing the cell-by-cell variability of the biological response discovered through use of this software.

# How to install:
1. Download the [contents of this repository](https://docs.github.com/en/get-started/start-your-journey/downloading-files-from-github#downloading-a-repositorys-files) and save the confocal_image_quantification folder including its contents on your machine.
2. Download [Ilastik](https://www.ilastik.org/)
3. Download pre-trained Ilastik segmentation project files for each image channel [here](https://doi.org/10.5281/zenodo.14537397) (you need the *.ilp files contained in that repository) or create your own. By default, the software expects to find these project files with the correct file-names in the same location as the .bat file used to run the software (details - see further on).
4. Download Python and the required dependencies (see ‘Dependencies’ section)

# Image Data Preparation and Software Usage

## 1. Prepare Your Image Data
- Separate the channels of interest into individual `.tif` files. Use [FIJI](https://imagej.net/software/fiji/downloads) with the [Olympus plugin](https://imagej.net/formats/olympus) to split the channels from Olympus `.oir` files.
- The required channels are:
  - **DNA signal ('blue')**
  - **Events of interest (e.g., endosomal escape) with general cell staining ('green')**
  - **Nuclei signal ('nuclei')**
- Save these `.tif` files in the same folder with the following names:
  - '(your_filename)**_blue.tif'**, '(your_filename)**_green.tif'**, '(your_filename)**_nuclei.tif'** (for 2D images)
  - '(your_filename)**_blu3.tif'**, '(your_filename)**_gr33n.tif'**, '(your_filename)**_nucl3i.tif'** (for Z-stacks)
- When analyzing a time-series of z-stacks (3D images), adjust the z_dimension_scale_factor in GUI_settings.ini (open in any text editor) to match the image acquisition parameters. The scale factor accounts for much smaller the resolution is in the z-axis as compared to the x and y. A factor of 0.1 would for example mean that the x/y image resolution is 10 times larger than the z-axis resolution.

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
