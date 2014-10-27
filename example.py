import transphere_python as tp

# will read ing example_input
example = tp.Transphere(modelfile='example_input_iras2a.inp') 
# at this point you can in IPython type "example." and 
# press "tab" to get various attributes that you can plot
# e.g. example.freq vs example.stellar_spec 
# for frequency vs. stellar spectrum. (BB-curve)
# or example.interp_opacities for the opacities

#then we write the input to files in the folder
example.write_input()

# annnd then we run it!
# you need transphere in your path for this to work.
example.run()
