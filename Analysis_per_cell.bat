::MIT License

::Copyright (c) 2024 Samuel Waldner s.waldner@unibas.ch

::Permission is hereby granted, free of charge, to any person obtaining a copy
::of this software and associated documentation files (the "Software"), to deal
::in the Software without restriction, including without limitation the rights
::to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
::copies of the Software, and to permit persons to whom the Software is
::furnished to do so, subject to the following conditions:

::The above copyright notice and this permission notice shall be included in all
::copies or substantial portions of the Software.

::THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
::IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
::FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
::AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
::LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
::OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
::SOFTWARE.

@echo off
setlocal enabledelayedexpansion

set "mainDirectory=%~dp0"
set IlastikProjectBlue=%mainDirectory%\pixel_classificaiton_blue.ilp
set IlastikProjectGreen=%mainDirectory%\pixel_classificaiton_green.ilp
set IlastikProjectGreenCells=%mainDirectory%\pixel_classificaiton_green_cells.ilp
set IlastikProjectNuclei=%mainDirectory%\pixel_classificaiton_nuclei.ilp
set IlastikExecutable=C:\Program Files\ilastik-1.4.0-gpu\ilastik.exe
set IlastikProjectBlue3d=%mainDirectory%\pixel_classificaiton_blue_3d.ilp
set IlastikProjectGreen3d=%mainDirectory%\pixel_classificaiton_green_3d.ilp
set IlastikProjectGreenCells3d=%mainDirectory%\pixel_classificaiton_green_cells_3d.ilp
set IlastikProjectNuclei3d=%mainDirectory%\pixel_classificaiton_nuclei_3d.ilp

for /r "%mainDirectory%" %%i in (.) do (
    set "subDirectory=%%i"
    set "subDirectory=!subDirectory:~0,-1!"
    if not "!subDirectory!"=="%mainDirectory%" (
        for %%j in ("!subDirectory!\*.tif") do (
            set "tiffFile=%%~nj"
            set "segmentedFile=!subDirectory!!tiffFile!_segmented.tiff"
            echo !tiffFile!
            echo !subDirectory!
            echo !segmentedFile!
            if not exist "!segmentedFile!" (
                echo "!tiffFile!" | find /i "blue" > nul && (
                    echo Segmenting blue spots using Ilastik for: "!subDirectory!" and "!tiffFile!"
                    "%IlastikExecutable%" --headless --output_format="multipage tiff" --export_source="Simple Segmentation" --output_filename_format={dataset_dir}\{nickname}_segmented.tif --project="%IlastikProjectBlue%" "!subDirectory!!tiffFile!.tif"
                ) || (
                    echo "!tiffFile!" | find /i "green" > nul && (
                        echo Segmenting green spots using Ilastik for: "!subDirectory!" and "!tiffFile!"
                        "%IlastikExecutable%" --headless --output_format="multipage tiff" --export_source="Simple Segmentation" --output_filename_format={dataset_dir}\{nickname}_segmented.tif --project="%IlastikProjectGreen%" "!subDirectory!!tiffFile!.tif"
                        echo Segmenting cells using Ilastik for: "!subDirectory!" and "!tiffFile!"
                        "%IlastikExecutable%" --headless --output_format="multipage tiff" --export_source="Simple Segmentation" --output_filename_format={dataset_dir}\{nickname}_cells_segmented.tif --project="%IlastikProjectGreenCells%" "!subDirectory!!tiffFile!.tif"
                    ) || (
                        echo "!tiffFile!" | find /i "nuclei" > nul && (
                            echo Segmenting nuclei using Ilastik for: "!subDirectory!" and "!tiffFile!"
                            "%IlastikExecutable%" --headless --output_format="multipage tiff" --export_source="Simple Segmentation" --output_filename_format={dataset_dir}\{nickname}_segmented.tif --project="%IlastikProjectNuclei%" "!subDirectory!!tiffFile!.tif"
                        ) || (
                            echo "!tiffFile!" | find /i "blu3" > nul && (
                                echo Segmenting 3D DNA spots using Ilastik for: "!subDirectory!" and "!tiffFile!"
                                "%IlastikExecutable%" --headless --output_format="multipage tiff" --export_source="Simple Segmentation" --output_filename_format={dataset_dir}\{nickname}_segmented.tif --project="%IlastikProjectBlue3d%" "!subDirectory!!tiffFile!.tif"   
                            ) || (
                                echo "!tiffFile!" | find /i "gr33n" > nul && (
                                    echo Segmenting 3D green spots using Ilastik for: "!subDirectory!" and "!tiffFile!"
                                    "%IlastikExecutable%" --headless --output_format="multipage tiff" --export_source="Simple Segmentation" --output_filename_format={dataset_dir}\{nickname}_segmented.tif --project="%IlastikProjectGreen3d%" "!subDirectory!!tiffFile!.tif"
                                    echo Segmenting 3D cells using Ilastik for: "!subDirectory!" and "!tiffFile!"
                                    "%IlastikExecutable%" --headless --output_format="multipage tiff" --export_source="Simple Segmentation" --output_filename_format={dataset_dir}\{nickname}_cells_segmented.tif --project="%IlastikProjectGreenCells3d%" "!subDirectory!!tiffFile!.tif"
                                ) || (
                                    echo "!tiffFile!" | find /i "nucl3i" > nul && (
                                        echo Segmenting 3D nuclei using Ilastik for: "!subDirectory!" and "!tiffFile!"
                                        "%IlastikExecutable%" --headless --output_format="multipage tiff" --export_source="Simple Segmentation" --output_filename_format={dataset_dir}\{nickname}_segmented.tif --project="%IlastikProjectNuclei3d%" "!subDirectory!!tiffFile!.tif"
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )
)

python "Colocalization_quantification_per_cell.py"

pause
endlocal
