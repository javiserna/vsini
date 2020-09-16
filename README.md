# Vsini by Fourier's Method

© Javier Serna, Jesús Hernández.

This is a semi-automatic tool designed to estimate the rotational velocity (vsini), using spectral lines of high resolution spectra. The code were developed in Python 3.7 and the GUI in Qt5.

[Hey! look at me](https://docs.google.com/presentation/d/1Cp8NaBN0EEg1mPtAIRAs3g8YdRQxtVRPUWPVnEFa0lI/edit?usp=sharing).

### New Features!

  - Gausian Fitting is done for the spectral line of interest (Automatic selection of the center and width of the line).

  - Automatic saving at the event file "Fourier.out".  The output file contains: (1) The name of the input file, (2) date, (3) time, (4) line midwave, (5) line width, (6) Spectral resolution, (7) Limb-Darkening coefficient, (8) vsini in km/s, (9) vsini uncertainty in km/s.
  - Plot the histrogram of vsini measures.

### Requirements and Dependencies

This tool is compatible with Python>2.7 on Linux, Mac y Windows. Some packages are necessary for their proper operation.

* Pandas
* Matplotlib.
* Tk
* Scipy
* PyQt5
* LMFIT

In case of not having any dependency, we suggest to install it in this way:

```sh
$ sudo apt install python-pandas
$ sudo apt install python-matplotlib
$ sudo apt install python-tk
$ sudo apt install python-scipy
$ sudo apt install python-pyqt5
$ sudo apt install python-lmfit
```

If you have an environment for Python> 2.7 (e.g, anaconda). The installation could be done in this way:

```
~(env)$ pip install pandas
~(env)$ pip install matplotlib
~(env)$ pip install tk
~(env)$ pip install scipy
~(env)$ pip install pyqt5
~(env)$ pip install lmfit
```

### How to use it?

Execute:

```sh
$ python vsini.py
```

The graphical interface of the program will appears

<img src="/home/javier/Desktop/fourier.png" alt="Interfaz Gráfica" style="zoom:80%;" />

By default, the spectral resolution (R=22500) and limb darkening coefficient (0.6) are optimally set up to work with APOGEE spectra. These can be changed by as the user requires. 

#### Step 1

Please press the "Load" button to get the spectra. CSV format for the input file is required, two columns of information, "col1": Wavelength and "col2": Flux. (e.g, load our test spectra "spec_to_test.csv").

#### Step 2 (Spectral Line Selection)

Press "Plot Spectrum" button to plot the input spectra. 

Use the “o” key and left click to zoom in the line of interest, or press "zoom to rectangle" icon.



<img src="/home/javier/Pictures/Screenshot from 2020-08-06 18-30-17.png" style="zoom:80%;" />

<img src="/home/javier/Pictures/Screenshot from 2020-08-06 18-30-34.png" style="zoom: 50%;" />





On each side of the line, press double click (A double click display a vertical red line). When it is done, the program performs a Gaussian fit to the spectra within the defined window. It computes the midwave of the line and width as 3 sigma level. These parameters are updated to the text boxes  "Line Center" and "Line Width" into the main window, once the user press "q" key.

<img src="/home/javier/Desktop/lineselection.gif" style="zoom: 50%;" />

<img src="/home/javier/Pictures/Screenshot from 2020-08-06 18-39-58.png" style="zoom:50%;" />

#### Step 3

Once the spectral line is selected, a visual inspection is necessary to check any issue  (e.g, mixed lines, bad normalization, bad pixels, telluric lines, etc.). Press "Plot Line" button, to plot the spectral line that will be use to measure vsini. 

Once it is certain that the line does not present any problem. Lets continue to the next step (look at recommendations).

#### Last Step

By pressing the "Run" button, vsini estimation is made. At the end, the tool plots a histogram with all random estimations of vsini. by pressing "q" key, vsini and their uncertainty is loaded at the main window. Finally, all the information is saved at the event file "Fourier.out".



------

© Instituto de Astronomía / UNAM (Ensenada, B.C, Mexico)

15 September 2020
