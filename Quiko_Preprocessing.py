
import os
import librosa

import numpy as np
from math import pi
import matplotlib.pyplot as plt

from scipy import signal
from scipy.signal import butter, lfilter, freqz

from scipy.stats import rv_continuous

from QuiKo_Backend import *
from Quiko_Analysis import * 
from Quiko_Feature_extraction import *

class sin_prob_dist(rv_continuous):
    """ This code is based on the example code from PennyLane """
    def _pdf(self, theta):
        # The 0.5 is so that the distribution is normalized
        return 0.5 * np.sin(theta)

# Samples of theta should be drawn from between 0 and pi
sin_sampler = sin_prob_dist(a=0, b=np.pi)

def random_dbmatrix( num_bands, qubits ):
    """
    This is to be called outside this library, mainly for the database 
    """
    random_matrix = []
    for band in range( num_bands ):
        #for index, subdiv in enumerate( range(2**qubits) ):
        phi, omega = np.random.uniform(0, 2*np.pi,size=2) # Sample phi and omega as normal
        theta = sin_sampler.rvs(size=1) # Sample theta from our new distribution
        angles = [round(theta[0],3), round(phi,3), round(omega,3)]
        
        random_matrix.append( angles )
        #random_matrix.append( random_matrix_b )
    return random_matrix

## FILTER BANK For all audio, This code was based off code found on StackExchange...
def normalized(a, axis=-1, order=2):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)
## LOWPASS FILTER---------------------------------------------

def butter_lowpass( cutoff, fs, order=5 ):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5 ):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

## BANDPASS FILTER---------------------------------------------
def butter_bandpass( lowcut, highcut, fs, order=5 ):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter( data, lowcut, highcut, fs, order=5 ):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

## HIGHPASS FILTER---------------------------------------------
def butter_highpass( cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter( data, cutoff, fs, order=5 ):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y

## Utility Functions--------------------------------------------------------------------
def ifnan( harm_est, spec_cent, perc_est ):
    if np.isnan(harm_est) == True:
        harm_est = 2*np.pi
    if np.isnan(spec_cent) == True:
        spec_cent = 2*np.pi
    if np.isnan(perc_est) == True:
        perc_est = 2*np.pi
        
    return harm_est, spec_cent, perc_est

def ifinf( harm_est, spec_cent, perc_est ):
    if np.isinf(harm_est) == True:
        harm_est = np.pi
    if np.isinf(spec_cent) == True:
        spec_cent = np.pi
    if np.isinf(perc_est) == True:
        perc_est = np.pi
        
    return harm_est, spec_cent, perc_est

def parse_featuredict( database: list ) -> list:
    """
    This is mostly ued for initializing the database...
    This is currently harcoded for the 3 qubit subband case...This needs to be changed as
    the subbands need to incr with the subdivision qubit register 
    """
    feature_matrix = []; track_names = []
    for track_data in database:
        bundle = []; filenames = []
        track_title = [*track_data][0]
        track_names.append( track_title )
        subbands = []
        for bands in track_data[track_title]:
            f_set = []
            f_set.append( track_data[track_title][bands]['harm_est'] )
            f_set.append( track_data[track_title][bands]['spec_cent'] )
            f_set.append( track_data[track_title][bands]['perc_est'] )

            subbands.append( f_set )
        feature_matrix.append( subbands )

    return track_names, feature_matrix
##-------------------------------------------------------------------------------------
def execute_quantumcircuit( quantumcircuit, dev='aer_simulator'): #'ibmq_16_melbourne' ):
    """
    dev: this specifies which device to send the circuit to, or use simulator
         --> dev='aer_simulator'
    """
    if dev == 'aer_simulator':
        aer_sim = Aer.get_backend('aer_simulator')
        transpiled_quiko_circuit = transpile(quantumcircuit, aer_sim)
        qobj = assemble(transpiled_quiko_circuit)
        results = aer_sim.run(qobj).result()
        distribution = results.get_counts()
        #plot_histogram(distribution); plt.show()

    else:
        provider = IBMQ.get_provider(hub='ibm-q-academic', group='stanford', project='q-comp-music-gen')
        device = backend = provider.get_backend( dev ) 
        
        # Run our circuit on the least busy backend. Monitor the execution of the job in the queue
        from qiskit.tools.monitor import job_monitor

        shots = 1024
        transpiled_circuit = transpile( quantumcircuit, backend, optimization_level=3 )
        job = backend.run( transpiled_circuit )
        job_monitor( job, interval=2 )
        
        results = job.result()
        distribution = results.get_counts()
        #plot_histogram(distribution)
    
    return distribution
##-------------------------------------------------------------------------------------
def Init_Database( num_bands: int, band_mapping: dict, directory: str, backend, 
                   all_in_one=False, random_anal=False ):
    """
    This should generate the database automatically...The input should be a list of dicts
    structure: { filename: { low: {features}, mid: {features}, high{features} } }
    feature_dict2 = []
    This is if the audio sample database need to be initialized
    """
    features = sample_database( num_bands, band_mapping, directory, random_anal=random_anal )
    track_names, feature_matrix = parse_featuredict( features ) ## I do not think I need parsefunction
    #feature_dict2.append( feature_dict[0] ); print( feature_dict2 )
    #print('hi It is calling it again here if I appear...')
    
    #print( track_names )
    #print( feature_matrix )
    bundle = []; filenames = []
    for index, track in enumerate( track_names ):
        qkc = QuikoCircuit( num_bands, 0, 0, feature_matrix[index], Encoding_method=0, big_circ=True )
        quiko_circ = qkc.Quantum_Circuit( 0,[] )
        bundle.append(quiko_circ)

    if all_in_one == False:
        q_state = execute_quantumcircuit( bundle, dev=backend )
        for index, result in enumerate(q_state):
            q_state[index] = zero_fill( q_state[index], num_bands ) ## result will be the bundle array

        return bundle, q_state, track_names
    else:
        dataB_dict = {}
        for index, entry in enumerate( track_names ):
            dataB_dict[entry] = bundle[index] #and  I just Have to reference it now
        # I could just return it in the function...
        return dataB_dict # This is if I am preparing to compare in one circuit

##------------------------------------------------------------------------------------    
def scrap_samples( folder ):
    """
    This function is OBSELETE but keep just in case I need it again...
    """
    Drum_Samp = {} 
    for file in os.listdir(folder):
        if file.endswith( ".wav" ):
            path = './{}/{}'.format(folder, file)
            drum, fs = librosa.load(path)
            Drum_Samp[file] = [drum, fs]

    return Drum_Samp

## Note: I might need to organize into classes with a self attribute here...-----------
def samp_subband( num_bands: int, folder: str, band_mapping: dict ) -> dict:
    """This is just parsing out the subband array"""
    Drum_Samp = {}
    for file in os.listdir( folder ):
        if file.endswith( ".wav" ):
            path = './{}/{}'.format(folder, file)
            drum, fs = librosa.load(path, sr=44100)
            subband = subband_gen( drum, fs, num_bands, band_mapping )
            Drum_Samp[file] = [subband, fs]

    return Drum_Samp

def subband_gen( input_meas, fs, num_bands, band_mapping ):
    """
    generate the desired subbands and corresponding spinal cord register size
    """
    ## Keep ing mind I should have an error check, the number of bands needs to correspond
    ## to the number of band mappings in a dictionary { b0: , b1: , ..., bn: }
    ## Make sure to have an error check between bands and time beat resolution

    ## Note: That this Function can also be used for the input track
    subbands = []
    for band in range( num_bands ):
        if band == 0:
            #print('0')
            freq_low = band_mapping[band]
            prev_freq = freq_low

            input_meas_low = butter_lowpass_filter( input_meas, freq_low, fs )
            input_meas_low = np.nan_to_num(input_meas_low)
            input_meas_lowN = input_meas_low / np.sqrt(np.sum(input_meas_low**2))
            #print(input_meas_lowN)
            
            subbands.append( input_meas_low )

        elif ( band < num_bands-1 ) & ( band > 0 ):
            #print('check')
            freq_upband = band_mapping[band]
            input_meas_mid = butter_bandpass_filter( input_meas, prev_freq, freq_upband, fs )
            input_meas_mid = np.nan_to_num(input_meas_mid)
            input_meas_midN = input_meas_mid / np.sqrt(np.sum(input_meas_mid**2))
            
            prev_freq = freq_upband
            subbands.append( input_meas_midN )

        else:
            input_meas_high = butter_highpass_filter( input_meas, prev_freq, fs )
            input_meas_high = np.nan_to_num(input_meas_high)
            input_meas_highN = input_meas_high / np.sqrt(np.sum(input_meas_high**2))
            
            subbands.append( input_meas_highN )

    return subbands

def sample_database( num_bands: int, band_mapping: dict, directory: str, random_anal: bool ) -> list:
    """
    This function generates the feature matrix for each band in the tracks for the
    database
    """
    Drums = samp_subband( num_bands, './dir/{}'.format(directory), band_mapping ) ## this is a dictionary

    features = []
    labels = Drums.keys()
    for label in labels:
        fs = Drums[label][1]
        subbands = Drums[label][0]

        sources = p_proc( subbands, fs)

        ## I need to now fix the feature_map function to go beyond 3 bands
        audio_feat = feature_map( label, sources, fs, random_anal ) #-> this should be a dictionary
        features.append(audio_feat)

    return features

## Feature Mapping for Database -------------------------------------------------------------
def feature_map( label, sources, fs, rand_anal ) -> list:
    feature_set  = {}; feature_setB = {}; feature_map = {}
    #rand_m = random_dbmatrix( len(sources), 3 ) ## Hardcode 3 qubits for now...need unhardcode it here NOW!
    rand_m = random_dbmatrix( len(sources), len(sources) )
    #print(sources)
    #print( ' num_bands: ', len(sources) )
    for band_num, band in enumerate( sources ):
        ## This needs to be revised...
        if rand_anal == False:
            band_data = sources[band]
            harm_est0  = harm(band_data[0], fs)
            spec_cent0 = spectral_centroid(band_data[0], fs) 
            perc_est0  = perc(band_data[1], fs)

            harm_est = 360 / harm_est0
            spec_cent = 360 / spec_cent0
            perc_est = perc_est0 / 360 #(2*np.pi) # This needs to be scaled down so the phase kickback is not so large
            
            harm_est, spec_cent, perc_est = ifnan( harm_est, spec_cent, perc_est )
            harm_est, spec_cent, perc_est = ifinf( harm_est, spec_cent, perc_est )

        else:
            harm_est  = rand_m[band_num][0]
            spec_cent = rand_m[band_num][1]
            perc_est  = rand_m[band_num][2]
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        feature_set['harm_est'] = round( harm_est, 3 )
        feature_set['spec_cent']=  round( spec_cent, 3 )
        feature_set['perc_est'] = round( perc_est, 3 )
        #print( 'dict check: ', feature_set)
        feature_setB[band_num] = feature_set
        #print( 'master_dict: ', feature_setB )

        feature_map[label] = feature_setB
        feature_set = {}
        
    return feature_map

def feature_dict( features_dict: list ) -> list:
    """
    Just to make sure that the formatting for the feature_matrix is correspoding with the backend...
    I still need to use it later for dictionary with the feature_matrix being the value, filename being the key
    """
    feature_matrix = []
    for band, audio_features in enumerate( features_dict ):
        band_seq = []
        for subdiv_val, subdiv_key in enumerate(audio_features.keys()):
            subdiv_feat = []
            subdiv_feat.append( audio_features[subdiv_key]['harm_est']  )
            subdiv_feat.append( audio_features[subdiv_key]['spec_cent'] )
            subdiv_feat.append( audio_features[subdiv_key]['perc_est']  )

            band_seq.append( subdiv_feat )
        feature_matrix.append( band_seq )

    return feature_matrix