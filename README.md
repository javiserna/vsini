# Vsini by Fourier's Method

© Javier Serna, Jesús Hernández and ARYSO group.

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fjaviserna%2Fvsini&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=Visits&edge_flat=false)](https://hits.seeyoufarm.com)

Please cite:
[![DOI:10.3847/1538-4357/AC300A](https://badgen.net/badge/DOI/10.3847/1538-4357/AC300A/blue/blue)](https://doi.org/10.3847/1538-4357/AC300A)

This is a semi-automatic tool designed to estimate the rotational velocity (vsini), using spectral lines of high resolution spectra. The code were developed in Python 3.7 and the GUI in Qt 5.

[Overview](https://docs.google.com/presentation/d/1Cp8NaBN0EEg1mPtAIRAs3g8YdRQxtVRPUWPVnEFa0lI/edit?usp=sharing).

### TO DO
- Add a new column with the relative error in to Fourier.out
- Pop up notifying error in line fitting. Suggest to choose other line or report as a not detection lines.

### Features!

  - Gausian Fitting is done for the spectral line of interest (Automatic selection of the center and width of the line).

  - Automatic saving at the event file "Fourier.out".  The output file contains: (1) The name of the input file, (2) date, (3) time, (4) line midwave, (5) line width, (6) Spectral resolution, (7) Limb-Darkening coefficient, (8) vsini in km/s, (9) vsini uncertainty in km/s.
  - Plot the histrogram of vsini measures.

### Requirements and Dependencies

This tool is compatible with Python>2.7 on Linux, Mac y Windows. Some packages are necessary for their proper operation.

* git
* Pandas
* Matplotlib.
* Tk
* Scipy
* PyQt5
* LMFIT

In case of not having any dependency, we suggest to install it in this way:

```zsh
$ sudo apt install git
$ sudo apt install python-pandas
$ sudo apt install python-matplotlib
$ sudo apt install python-tk
$ sudo apt install python-scipy
$ sudo apt install python-pyqt5
$ sudo apt install python-lmfit
```

If you have an environment for Python> 2.7 (e.g, anaconda). The installation could be done in this way:

```zsh
~(env)$ pip install pandas
~(env)$ pip install matplotlib
~(env)$ pip install tk
~(env)$ pip install scipy
~(env)$ pip install pyqt5
~(env)$ pip install lmfit
```

### How to use it?

First of all, we are going to create a folder where we will download the package:

```zsh
$ mkdir Fourier
$ cd Fourier
$ git clone https://github.com/javiserna/vsini
$ cd vsini
```

and Execute:

```zsh
$ python vsini.py
```

The graphical interface of the program will appears

<img src="https://raw.githubusercontent.com/javiserna/vsini/master/Images/fourier.png?token=ADW2GZ3M4Z46V4WHSQUX63K7MGPVA" alt="Interfaz Gráfica" style="zoom:80%;" />

By default, the spectral resolution (R=22500) and limb darkening coefficient (0.6). These can be changed by as the user requires. 

#### Step 1

Please press the "Load" button to get the spectra. CSV format for the input file is required, two columns of information, "col1": Wavelength and "col2": Flux. (e.g, load our test spectra "spec_to_test.csv").

#### Step 2 (Spectral Line Selection)

Press "Plot Spectrum" button to plot the input spectra. 

Use the “o” key and left click to zoom in the line of interest, or press "zoom to rectangle" icon.



<img src="https://raw.githubusercontent.com/javiserna/vsini/master/Images/Screenshot%20from%202020-08-06%2018-30-17.png?token=ADW2GZ2NCEBW2W7C3BG4VPK7MGPZ2" style="zoom:80%;" />

<img src="https://raw.githubusercontent.com/javiserna/vsini/master/Images/Screenshot%20from%202020-08-06%2018-30-34.png?token=ADW2GZY57FEETFLBAB2BKWS7MGP4I" style="zoom: 50%;" />





On each side of the line, press double click (A double click display a vertical red line). When it is done, the program performs a Gaussian fit to the spectra within the defined window. It computes the midwave of the line and width as 3 sigma level. These parameters are updated to the text boxes  "Line Center" and "Line Width" into the main window, once the user press "q" key.

<img src="https://raw.githubusercontent.com/javiserna/vsini/master/Images/lineselection.gif?token=ADW2GZ4KOZNTYOGAQ3MBIYC7MGP52" style="zoom: 50%;" />

#### Step 3

Once the spectral line is selected, a visual inspection is necessary to check any issue  (e.g, mixed lines, bad normalization, bad pixels, telluric lines, etc.). Press "Plot Line" button, to plot the spectral line that will be use to measure vsini. 

Once it is certain that the line does not present any problem. Lets continue to the next step (look at recommendations).

#### Last Step

By pressing the "Run" button, vsini estimation is made. At the end, the tool plots a histogram with all random estimations of vsini. by pressing "q" key, vsini and their uncertainty is loaded at the main window. Finally, all the information is saved at the event file "Fourier.out".

### Checking the Final Results

Once we have finished our measurements, then close the window and execute: 

```zsh
$ python results.py
```
This script will organize the historical measurements located in the file "Fourier.out" For the set of spectra imported and the set of lines processed by the user. Finally, the user will have a file called "results.csv" with the spectra name, vsini, and vsini uncertainty.

------

© Instituto de Astronomía / UNAM (Ensenada, B.C, Mexico)

15 September 2020
