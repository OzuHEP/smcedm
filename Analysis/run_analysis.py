#!usr/bin/python

import os, sys

import configparser

config_file_name = sys.argv[1]

config = configparser.ConfigParser()
config.read(str(config_file_name))

infile_name   = str(config['io']['infile_name'])
outfile_name  = str(config['io']['outfile_name'])
n_btag        = str(config['cuts']['n_btag'])
n_light_jet   = str(config['cuts']['n_light_jet'])
MET_threshold = str(config['cuts']['MET_threshold'])

print("analysis is running with the following parameters:")
print( "input file_name: " + infile_name + "\n" + "output file_name: " + outfile_name + "\n" +  "number of required btag jets: " + n_btag + "\n" + "number of required light jets: " + n_light_jet + "\n" + "MET threshold: " + MET_threshold + " GeV")

os.system( "./SMCEDM_Analysis " + infile_name + " " + outfile_name + " " + n_btag + " " + n_light_jet + " " + MET_threshold)
os.wait()