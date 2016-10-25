import pypie
import matplotlib.pyplot as plt
import numpy as np

from periodictable import formula

# Set Plot Stlye
plt.style.use('ggplot')

# Load Data and Create pypie.Pie object
data = 'sample_data' 
pie = pypie.Pie(data)

energy = 13.1

# Calibrate Mass Spectra
pie.ms_cursor()
m1, m2 = formula('He').mass, formula('O2').mass
t1, t2 =  825, 3201 
pie.ms_calibrate( m1, m2, t1, t2)  # terms[file_key] format is ((m1, m2), (t1, t2))

# Plot and save mass spectra
pie.ms_plot() 
save_path = 'sample_mass_spectrum.txt'
pie.ms_save(path=save_path)


# Slice multiple PIEs with given formulae
formulae = [ 'CH4','C2H2','C2H4','C3H6','C4H6','C5H10','C6H6','C7H8','C7H8','C8H10']
masses = [ formula(form).mass for form in formulae ]
for f, m in zip(formulae, masses):
        pie.pie_slice(m, 0.5, label=f)

# Show PIE info
pie.pie_info()

# Optional Signal Corrections
pie.current_plot()
pie.pie_current_correction()
'''
pie.pie_background_correction(background)
pie.pie_normalize()
'''

# Plot and Save PIEs
pie.pie_plot() 
pie.pie_save()
