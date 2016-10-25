# pypie
pypie is a 3D data processing library for rapid calibration, reduction, and presentation of mass spectrometric data series.

## Motivation

Built on python, pypie uses numpy and matplotlib to calibrate, reduce and plot 3D mass spectrometric data, especially those data recorded as a function of **P**hoto**I**onization **E**fficiency, hence py**pie**. It was built for use with the data generated at the [Chemical Dynamics Beamline](http://www.chemicaldynamics.lbl.gov/) 9.0.2, Molecular Beams Machine, LabView data acquisition software. It treats data where m/z vs intensity are acquired as function of synchrotron energy. However, it may be used or adapted for any mass spectrometric data series with multiple mass spectra over a given function, such as mass spectrum vs time.

## Dependencies
* python 2.7+
* numpy
* matplotlib
* periodictable

## Usage

pypie, `pypie.py`, may be run interactively from the python command line or invoked within a script for massive data throughput. `Pie` may also be inherited with modification of the `load` method for use with your own data types. 

In this example we will make use of [sample_data](http://) for demonstration of it's use. The sample_data is a series of mass spectra over a given energy range measured from calibrant gas containing helium and trace hydrocarbons, each at 100ppm. These hydrocarbons are: hydrogen, methane, acetylene, ethene, propene, 1,3-butadiene, cyclopentane, benzene, toluene, para-xylene. CHeck out the quickstart section for an overview of the commands to generate, plot, and save PIEs. A more thorough description of each method follows.

### Quickstart

```python
import pypie

from periodictable import formula

# Create pypie object from data
data = '/path/to/sample_data'
pie = pypie.Pie(data)

# Mass calibration from data points
m1, m2 = formula('He').mass, formula('O2').mass
t1, t2 = 825, 3201
pie.ms_calibrate(m1, m2, t1, t2)

# Define formulae for slicing
formulae_for_slicing = [
    'H2',         # hydrogen
    'CH4',        # methane
    'C2H2',       # acetylene
    'C2H4',       # ethene
    'C3H6',       # propene
    'C4H6',       # 1,3-butadiene
    'C5H10',      # cyclopentane
    'C6H6',       # benzene
    'C7H8',       # toluene 
    'C8H10',      # para-xylene
    ]

# Define masses from formulae
masses = [ formula(form).mass for form in formulae ]

# Slice spectrum at each mass and a given width to generate PIEs
for formula, mass in zip(formulae, masses):
    pie.pie_slice(mass=mass, width=0.5, label=formula)

# Plot slices on mass spectrum
pie.pie_show_slices()

# Plot current, current correction, and normalization (optional)
pie.current_plot()
pie.pie_current([])
pie.pie_norm([])

# Plot and Save PIEs
save_path = '/path/to/save/PIEs.dat'
pie.pie_save([])
pie.pie_plot([])

```
### Step-by-step Description

#### Creating a Pie object

First we import pypie and create a Pie object.

```python
import pypie

data = '/path/to/sample_data'
pie = pypie.Pie(data)

```

This will create a `Pie` object containing data read from `sample_data`. `sample_data` consists of a first column of time and subsequent columns of ion counts for each ionzation energy and photon current. In this case we measure a mass spectrum from 8 to 13.55 electronvolts (eV) with a step size of 0.075eV to give 74 columns of mass spectra. These data are stored as `pie.energy`, `pie.current`, `pie.time`, and `pie.counts`. We can list the available energies with `print(pie.energy)`. 

#### Calibrate the Mass Spectrum

After loading the data we need to calibrate its time component to m/z. To achieve this we should first inspect a mass spectrum at an energy with some known peaks to determine what those peaks correspond to in time (or data point in this case). This may be achieved with the `Pie.ms_cursor` which takes a specific energy as an argument such as one you pronted above. 

```python
pie.ms_cursor(13.1)
```
This will use matplotlib to plot the mass spectrum at 13.05eV. You should be familiar with your data to be able to identify which peaks might correspond to which mass. You will need two peaks to calibrate the m/z for the spectrum. In this case I know I have Helium (He) 825, and Oxygen (O2) at 3201. You can zoom on a particular peak with the zoom tool to get this precise value, the cursor will track your mouse position and display it next to the cursor. This is demonstrated below for the helium peak

![zoom tool](https://)

With these two times, and masses and using `Pie.cal_mass` method we can calibrate our time data to m/z. Using the `periodictable` library `formula` function makes defining the masses significantly easier since you can simply define the formula and infer the mass.

```python
from periodictable import formula

m1, m2 = formula('He').mass, formula('O2').mass
t1, t2 = 825, 3201

pie.ms_calibrate(m1, m2, t1, t2)
```
You can check the success of your calibration by plotting a mass spectrum at a particular energy.

```python
pie.ms_plot(13.1)
```
If you're not happy, find new terms with with `pie.ms_cursor` and recalibrate with `pie.ms_calibrate`. Once you're happy we can start extracting data over the energy mass series.

#### Slicing the MS Data vs Energy

We can slice data over the energy series using the `Pie.pie_slice(center, width, label)` where `center` and `width` are in units of m/z (Da), and `label` is a string. Again, this the function `formula` we imported above comes in handy. Remember gas mixture contained the hydrocarbons listed above. Taking toluene for example we can use it's molecular formula C7H8 to specifty as it's mass with formula using `formula('C<sub>7</sub>H8').mass` to return `92.13842`. The `width` of a given peak can be inspected with the `pie.ms_plot` method and the label can be a label of your choosing, or as I like to do, use the formula. Let's give it a go

```python
toluene = 'C7H8'
mass = formula(toluene).mass
width = 0.5
pie.pie_slice(mass, width, tolene)
```

Now we can visualize that slice with the `Pie.pie_show_slices` method


