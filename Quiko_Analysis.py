"""
Outline Analysis
0) Random values in the circuit (but I think this falls under expressibility...) -> Can be combined with (2)
1) Fidelity between audio and drum samples
2) expressibility between the different methods (to be gathered in one plot)

*** Classical analysis - Some similarity metric between the audio samples
"""
import csv
import numpy as np
from math import pi

import itertools
from itertools import combinations

import matplotlib.pyplot as plt
from scipy.stats import rv_continuous

from QuiKo_Backend import *
from Quiko_Preprocessing import *

from qiskit import IBMQ
from qiskit.providers.aer.noise import NoiseModel
from qiskit import QuantumCircuit, Aer, assemble, transpile, execute
from qiskit.visualization import plot_bloch_multivector, plot_histogram, array_to_latex

## Will need to import the analysis library from here this will connect to the results from
## running the Quantum Circuit

print('Loading analysis testing version...')
print('Loading IBM Account...')

subB = {2:[1080], # this is the 2-qubit case
        3:[920, 3150], # This is for the 3-qubit case
        4:[630, 1720,  4400], # This is for the 4-qubit case
        5:[510, 1270, 2700, 6400], # This is for the 5-qubit case
        6:[400, 920, 1720, 3150, 6400], # This is for the 6-qubit case
        8:[300, 630, 1080, 1720, 2700, 4400, 7700], # This is for the 8-qubit case
        12:[200, 400, 630, 920, 1270, 1720, 2320, 3150, 4400, 6400, 9500], # This is for the 12-qubit case
        25:[100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500]} # This is for the 25-qubit case

## Utility, Database & Execution Functions ---------------------------------------------
def execute_quantumcircuit_ANALY( quantumcircuit, inject_noise=False, dev='aer_simulator'): #'ibmq_16_melbourne' ):
    """
    dev: this specifies which device to send the circuit to, or use simulator
         --> dev='aer_simulator'
    """
    if dev == 'aer_simulator':

        if inject_noise == True:
            provider = IBMQ.load_account()
            backend = provider.get_backend('ibmq_belem')
            noise_model = NoiseModel.from_backend(backend)
            result = execute(quantumcircuit, Aer.get_backend('qasm_simulator'),
                 noise_model=noise_model).result()
            distribution = result.get_counts(0)
            #print('check: ', distribution)
        else:
            aer_sim = Aer.get_backend('aer_simulator')
            transpiled_quiko_circuit = transpile(quantumcircuit, aer_sim)
            qobj = assemble(transpiled_quiko_circuit)
            results = aer_sim.run(qobj).result()
            distribution = results.get_counts()

        #plot_histogram(distribution); plt.show()
    else:
        print('IN this execute function...')
        provider = IBMQ.get_provider(hub='ibm-q-academic', group='stanford', project='q-comp-music-gen')
        device = backend = provider.get_backend( dev ) 
        
        # Run our circuit on the least busy backend. Monitor the execution of the job in the queue
        from qiskit.tools.monitor import job_monitor
        print('all set up...')
        shots = 1024
        transpiled_circuit = transpile( quantumcircuit, backend, optimization_level=3 )
        print('backend_running...')
        job = backend.run( transpiled_circuit )
        print( 'job monitor is supposed to be on!')
        job_monitor( job, interval=2 )
        
        results = job.result()
        distribution = results.get_counts()
        #plot_histogram(distribution)
    return distribution

def database_retrieval( database, index, analysis=False ):
    """
    This function retrieves specific audio samples form the database...
    will be filled in when I start generating the database, This is iterating on the csv files generated
    """
    ## Read from csv file, save to variable then find the distr I want...more paring required here...
    with open("QuiKo_Drum_samp.csv") as f:
        for line in f:
            key, val = ( line.split(', ') )
            q_distr[key] = int( val ) ## Need to test this code snippet out before I execute everything
            
    return q_distr

def zero_fill( distr, qubit_count ):
    """
    This is going to be a utility function for filling in the states
    that the probabilities are zero (Already wrote one just kick it...)
    """
    keys = np.sort(list(distr.keys()))
    #print('keys: ',keys)
    #print(keys)
    targ_index = np.arange( 2**qubit_count )
    targ_key = [format(i, '0{}b'.format(qubit_count)) for i in targ_index]
    #print('targ_key: ',targ_key)
    diff = np.setdiff1d(targ_key, keys)
    #print('diff: ', diff)
    for i in diff:
        #print(i,type(i))ß
        distr[i] = 0
    #print('what is wrong: ',distr.keys(), distr)
    return distr
    #    if distr[k]
"""
The next 3 functions should go into Beat_Construction.py
"""
def timbre_per_beat( distribution, spine_qubits ):
    ## This might need to be relocated to beat construction module
    input_timbre = []
    for subdivision in range( 2**spine_qubits ):
        timbre_state = {}
        subd_bin = format(subdivision, "0{}b".format( spine_qubits ))
        
        for i in distribution.keys(): #list(distribution.keys()):
            #print( i )
            if subd_bin in i[spine_qubits:]:
                if i[:spine_qubits] in timbre_state:
                    timbre_state[i[:spine_qubits]] += distribution[i]
                else:
                    timbre_state[i[:spine_qubits]] = distribution[i]

        input_timbre.append( timbre_state )

    for index, div in enumerate( input_timbre ):
        input_timbre[index] = zero_fill(div, spine_qubits)

    return input_timbre

def PK_Sequence( distribution, spine_qubits ):
    """ 
    This function will show the effects of the PKSBE (PKBSE ONLY)

    Steps:
    1) Execute this function and get the distribution of subdivision codes
    2) Measure the equal superposition of the subdivisions
    3) find the fidelity between the two distributions
    """
    PK_seq = {}
    spine_codes = []
    
    for count, string in enumerate( distribution ):
        spine_codes.append( string[3:] )

    for s_code in np.unique(spine_codes):
        PK_seq[s_code] = 0
        for key in distribution.keys():
            if s_code in key:
                PK_seq[s_code] += distribution[key]
        
    return PK_seq


def SPP_Rhythm( distribution, spine_qubits, timbre_qubits ):
    """
    This is generate a superposition of rhythms...I might need a different analysis
    or adapted adapted expressibility to analyze this decoded beat array...

    --> This might need to be relocated to the beat construction module
    """
    ## 1) Do the same process in timbre_per_ but do it with respect to time
    ## 2) The difference here is that for the spine we don't want to see the timbre distr
    ##    We are trying to see the percentage that the timbre code appears on the subdivision set
    input_spine = []
    for timbre_code in range( 2**timbre_qubits ):
        spine_state = {}
        timb_bin = format( timbre_code, "0{}b".format(timbre+_qubits) )

        for i in distribution.keys():
            if timb_bin in i[:timbre_qubits]: ## The other way/the timbre register instead...
                if i[spine_qubits:] in spine_state:
                    spine_state[i[spine_qubits:]] += distribution[i] ## This part is a little different it should just be adding one
                else:
                    spine_state[i[spine_qubits:]] = distribution[i] # Ok this was kinda the same as the timbre_per_beat function
        input_spine.append( spine_state )

    for index, timbre in enumerate( input_spine ):
        input_spine[index] = zero_fill( timbre, timbre_qubits )
    
    return input_spine

""" End of the functions that need to be moved to Beat_Construction.py"""  
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

def bundle_qcircuits( qcirc_array: list ):
    """
    This is just a utility function to bundle all the circuits that I have for one job...
    I think I can just feed a list to the execute circuit function above...Keep this function here just in case
    """

    pass

##---------------------------------------------------------------------------- 
class fidelity_measure:
    """
    This class will be used to determine the audio and drums samples to be used
    from the database based on the quantum state of the input
    """
    def __init__( self, qubits, distr_1: dict, distr_2: dict,  init_database=True ):
        self.qubits = qubits
        self.distr_1 = distr_1
        self.distr_2 = distr_2
        
    def measured_qstate( self ):
        """
        returns the statevector and density matrix of the measured circuits
        This is mostly a parsing function... Keep in mind that i need implement the z fill for the actual thing
        """
        
        distr_1 = zero_fill( self.distr_1, self.qubits )
        distr_2 = zero_fill( self.distr_2, self.qubits )
        q_state1 = np.linspace(0, len(list(distr_1.keys())), num=len(list(distr_1.keys())))
        q_state2 = np.linspace(0, len(list(distr_2.keys())), num=len(list(distr_2.keys())))
        
       # print( 'distr_1', distr_1 )
       # print( 'distr_2: ', distr_2 )
        
        for index, state in enumerate( distr_1 ):
           #print('state: ', state)
            q_state1[index] = distr_1[state]
            
        for index, state in enumerate( distr_2 ):
            q_state2[index] = distr_2[state]
        
        q_state1 = q_state1 / np.max(q_state1)
        q_state2 = q_state2 / np.max(q_state2)
        return q_state1, q_state2
        
    def fidelity_meas( self ) -> float:
        """
        The fideltiy metric will be used to detemrine the closeness between
        the input and database quantum states for the timbre qubits
        
        I can also test this real quick with two state vectors that I know are close and distant
        """
        
        ## looking for the squared over lap of the two states
        q_state1, q_state2 = self.measured_qstate() ## Also these are just column vectors (kets)
        
        ro = q_state1 # Make sure they are of the type -> numpy arrays
        #print('q_state',q_state2)
        omega = q_state2 #.reshape( 2**self.qubits,1 )
        #fidelity = abs( ro * omega )**2
        fidelity = abs( np.inner(ro,omega) )**2
        #print('check: ', fidelity)
        return fidelity
    
##----------------------------------------------------------------------------    
class Expressibility:
    #def __init__( self, quantumcircuit, timbre_qubits, spine_qubits, qubits,
    #                    database_states: list ):#distr_1: dict, distr_2: dict ):
    def __init__( self, num_bands, timbre_qubits, spine_qubits, qubits, database_bundle): #distr1: list, database_states: list ):#distr_1: dict, distr_2: dict ):
        #self.qc = quantumcircuit
        
        self.qubits = qubits
        self.num_bands = num_bands
        self.spine_qubits = spine_qubits
        self.timbre_qubits = timbre_qubits
        
        self.dbb = database_bundle
        
        #self.distr_1sim = distr1; 
        #self.qstates = database_states #distr2 is the database states
        
    def Haar_meas( self, fidelities: dict ) -> dict:
        """
        This will be taken from then PennyLane example...
        I think the difference between that and the slides is how they treat the theta
        
        This is just an even distribution of fidelities depending on the number of unique fidelities
        uniq_fid -> is the number of fidelity measurements that the random_meas() produced
        """
        """
        size = 0
        hm = {}
        L = np.arange( 61 )
        for l in combinations(L,5):
            hm
        """
        haar_meas = {}; uniq_fid = len(fidelities)
        prob = 1 / uniq_fid
        
        for f in fidelities:
            haar_meas[f] = prob
        
        return haar_meas
  
    def random_meas( self, M_samp=200, encoding_method='pkbse', fid_meas=False, big_circ=False, dev='aer_simulator' ):
        """
        This will be taken from then PennyLane example...
        I think the difference between that and the slides is how they treat the theta
        
        (1) Ok the idea for this experiment is to generate random matrices and compare them t the database, Theoretically I should get 
            a flat distrubution of audio samples to be chosen from the database.
        
        (2) Samething as the expressibility measure but instead of the haar measure we are going to use the comparison between the input and
            output as the reference instead of the haar random measure. from there calaculate the expressibility in terms of the database
            
            *** I can proabably do this for timbre and the spinal cord qubits -> two different graphs
            *** This is going to be adapted so that the input state is random and is compared to the actual states of the 
                just look at the fidelities and see which fidelities appear more than others ( distrbution ) --> compare it to the haar random distribution
                
            *** Technically since the dataase samples are using generic u gates the this is the haar random measure
            
            What does this tell us -> What the distribution of exploring the audio sample space...the database -> Then do the process again with multiple dimensions (x,z,and y) basis larger circuit involved
        """
        """
        num_bands = 3
        band_mapping = {}
        
        for freq in range( len(subB[num_bands]) ):
            band_mapping[freq] = subB[num_bands][freq]
        band_mapping
        database_bundle = Init_Database( self.num_bands, band_mapping, 'Drums_Database', 'none', all_in_one=True, random_anal=True)
        """
        
        bundle = []
        """
        if fid_meas == True:
            M_samp = 1

        else:
            if 2**self.num_bands < M:
                M_samp = M - 2**self.num_bands
            else:
                M_samp = 2**self.num_bands
        """

        window_Rmeas = []; window_Hmeas = []
        # Generate the randomized matrix and collect M samples
        #for window in range( len(self.qstates) - 5 ): # This 5 value is the layer size!
        random_meas = {} 
        for M in range( M_samp ):
            random_matrix = []
            # Set up the feature matrix -> encoding matrix process
            for band in range( self.num_bands ):
                random_matrix_b = [] 
                for index, subdiv in enumerate( range(2**self.timbre_qubits) ):
                    theta, phi, omega = 2 * np.pi * np.random.uniform(size=3) # Sample phi and omega as normal
                    #theta = sin_sampler.rvs(size=1) # Sample theta from our new distribution
                    angles = [round(theta,3), round(phi,3), round(omega,3)]
                    
                    random_matrix_b.append( angles )
                random_matrix.append( random_matrix_b )
            
            # build and execute the circuit
            #x = Encoding_matrix( random_matrix, 8, 3, 3 )

            #print( len(random_matrix) )
            print( 'Sample: ', M )
            if encoding_method == 'pkbse':
                
                etype = 1
                x = Encoding_matrix( random_matrix, 2**self.num_bands, self.num_bands, 2**self.num_bands ) #Encoding_matrix( random_matrix, 8, 3, 3 )
                y = x.pkbse_encoding_matrix()
    
                if big_circ == False:
                    meas_m = np.arange( 2*self.num_bands ) #This is hard-coded for the case of num_bands = spine_qubits
                    cbits = self.timbre_qubits + self.spine_qubits
                    
                    qkc = QuikoCircuit( self.num_bands, self.spine_qubits, cbits, y, etype )
                    qc  = qkc.Quantum_Circuit( etype, meas_map=[meas_m,meas_m] )
                    """Right Here is where I'm going to have to call to larger circuit gen function """
                    bundle.append( qc )
                    #distr_1 = execute_quantumcircuit_ANALY( qc, dev) #inject_noise=True, dev=dev )
                    #distr_1 = timbre_per_beat( distr_1, self.spine_qubits )
                   
                else:
                    meas_m = np.arange( 2*self.num_bands ) #This is hard-coded for the case of num_bands = spine_qubits
                    cbits = self.timbre_qubits + self.spine_qubits
                    
                    qkc = QuikoCircuit( self.num_bands, self.spine_qubits, cbits, y, etype, big_circ=True )
                    qc  = qkc.Quantum_Circuit( etype, meas_map=[])#[meas_m,meas_m] )
                    
                    qc_full = Full_Quantum_Circuit( qc, self.dbb, 1 )
                    #print(' check: ', qc_full)
                    bundle.append( qc_full )

            elif encoding_method == 'static':
                etype = 2
                x = Encoding_matrix( random_matrix, 2**self.num_bands, self.num_bands, 2**self.num_bands )
                y = x.static_encoding_matrix()
                
                if big_circ == False:
                    meas_m = np.arange( 2*self.num_bands ) #This is hard-coded for the case of num_bands = spine_qubits
                    cbits = self.timbre_qubits + self.spine_qubits
                    
                    qkc = QuikoCircuit( self.num_bands, self.spine_qubits, cbits, y, etype )
                    qc  = qkc.Quantum_Circuit( etype, meas_map=[meas_m,meas_m] )
                    """Right Here is where I'm going to have to call to larger circuit gen function """
                    bundle.append( qc )
                    #distr_1 = execute_quantumcircuit_ANALY( qc, dev) #inject_noise=True, dev=dev )
                    #distr_1 = timbre_per_beat( distr_1, self.spine_qubits )
                    
                else:
                    meas_m = np.arange( 2*self.num_bands ) #This is hard-coded for the case of num_bands = spine_qubits
                    cbits = self.timbre_qubits + self.spine_qubits
                    
                    qkc = QuikoCircuit( self.num_bands, self.spine_qubits, cbits, y, etype, big_circ=True )
                    qc  = qkc.Quantum_Circuit( etype, meas_map=[])#[meas_m,meas_m] )
                    
                    qc_full = Full_Quantum_Circuit( qc, self.dbb, 1 )
                    #print(' check: ', qc_full)
                    bundle.append( qc_full )

            elif encoding_method == 'AP':
                """
                x = Encoding_matrix( random_dbmatrix( 3,3 ), 8, 3, 3)
            
                etype = 0
                y = x.APencoding_matrix()
                qkc = QuikoCircuit( 3, 3, 3, y, etype )
                qc  = qkc.Quantum_Circuit( etype, meas_map=[[0,1,2],[0,1,2]] )
                distr_1 = execute_quantumcircuit_ANALY( qc, dev='aer_simulator' )
                distr_1 = timbre_per_beat( distr_1, 3 )
                """
                distr_1 = self.distr_1sim
            else:
                print( 'Error: Invalid Encoding_Method, Pleas try again...' )
                return None
            
            ## Note we are ignoring the spinal qubit register in this case so we 
            # We only nned to measure the first 3 qubits (timbre register)
            
            #print('checking: ', distr_1)
            """right here i need to call a funcion that compares between the states in the database
            Will be done when the database is set up... """
            #distr_2 = database_retrieval( database, index ) ## Right here is where I need to revise
            ## qstates should be a radnomly generated database of a particular size
        print('Express working...')
        distr_1 = execute_quantumcircuit_ANALY( bundle, dev=dev) #inject_noise=True, dev=dev )

        import json
        with open('QuiKo_ANALY_3qb_larger', 'w') as fout:
            json.dump(distr_1, fout)

        return distr_1, self.dbb
    
    def Express_larger( self, distr_1, layer_start, layer_end):#layer_size ):
        ## This top section of this function needs to be a decoder function...
        track_titles = list( self.dbb.keys() )
        
        ## Generate the haar measure distrubution
        haar_meas = {}
        """
        L = np.arange( len(track_titles) )
        for l in combinations(L,5):
            haar_meas[str(list(l))] = 1
        """
        
        ## Generate the random_meas distribution
        result_database = []; rand_meas_arr = []
        for results in distr_1:
            keys = results.keys()
    
            detect = []
            spine  = []
            where  = []
            for i in keys:
                #print( len(i[:len(i)-3]) )
                detect.append(i[:len(i)-3] )
                spine.append( i[len(i)-3:] )
            uniq_subdiv = np.unique(spine)
            
            ##----------------------------
            markers = []
            for det in detect:
                init = -3
                samp_marker = []
                for index, compr in enumerate( range(int(len(det)/3)) ):
                    init += 3
                    frame = det[init:init+3]
                    #print('check: ', frame, index)
                    
                    cpr_flag = int(frame,2)
                    #print('check2: ', cpr_flag)
                    if cpr_flag == 0:
                        samp_marker.append(index)
                    #else:
                    #    sampe_marker.apppend(-1)
                markers.append(samp_marker)
            
            for j in uniq_subdiv:
                where.append( [i for i, e in enumerate(spine) if e == j] )
                
            ##----------------------------
            super_layer = {}
            for subdiv, loc in enumerate( where ):
                #print('for this subdiv: ', subdiv)
                layer = {}
                for shot in loc:
                    for sp in markers[shot]:
                        #print(sp)
                        
                        if sp in layer.keys():
                            layer[sp] += 1
                        else:
                            layer[sp] = 1
                super_layer[subdiv] = layer #We are actually talking about samples
                
            ##-----------------------------
            layer_struct = {}
            for subdiv in super_layer:
                layer_div = []
                distr_ = super_layer[subdiv]
                sort_orders = sorted(distr_.items(), key=lambda x: x[1], reverse=False)
                #print( 'beat: ', subdiv )
                for i in sort_orders[layer_start:layer_end]: # Here is where I specify the layer...
                    #print(i[0], i[1])
                    layer_div.append( i[0] )
                
                layer_struct[subdiv] = layer_div
            
            result_database.append( layer_struct ) # There should be 50 of these...
            
            random_distr = {}
            key_lst = []; key_label = -1
            
            for trial in result_database:
                for j in trial:
                    ptr_ = str(np.sort( trial[j] ))
                    if ptr_ in key_lst:#list(random_distr.keys()):
                        k = key_lst.index( ptr_ )
                        random_distr[k] += 1
                    else:
                        key_lst.append(ptr_)
                        k = key_lst.index( ptr_ )
                        random_distr[k] = 1
            #random_meas_arr.append( random_distr  )

        return random_distr, haar_meas

#express_p, pk_seqp = p.Expressibility_Measure( 'pkbse', 'ibmq_manhattan')
            
    def Express( self, distr_1, dev, fid_meas=False ):
        # layer = []
        #print('Express working...')
        #distr_1 = execute_quantumcircuit_ANALY( bundle, dev=dev) #inject_noise=True, dev=dev )

        """ 
        Right here I should have another function specifying the effect and behaviour of the phase kick
        back sequencing ...
        """
        pk_seq_per_samp = []
        window_Rmeas = []; window_Hmeas = []
        index = -1
        for window in range( len(self.qstates) ):
            random_meas = {}
            for samp in distr_1: # For each sample in the bundle...
                #print('new_distr: ', dist)

                for d1 in samp: # for each subdivision in the sample...
                    #print('new_distr_: ')

                    fidelity_array = np.linspace(0, len(self.qstates), num=len(self.qstates))
                    for track_id, state in enumerate( self.qstates ):
                        distr_pair = []
                        distr_pair.append( d1 )
                        distr_pair.append( state )
                        
                        F = fidelity_measure( self.qubits, distr_pair[0], distr_pair[1] )
                        fidelity_ = F.fidelity_meas()
                        
                        fidelity_array[track_id] = round( fidelity_, 5 )

                    layer = np.argsort( fidelity_array )
                    #print( layer )
                    ## Execute for the layer then the single audio sample case...
                    if fid_meas == True:
                        #layer = np.sort( fidelity_array )
                        return layer, fidelity_array

                    """ 
                    I can do some better analysis on the window layering.
                    Compare similarity between different layers... neighboring layers?
                    """
                    layer_dict ={}
                        
                    #index += 1
                    #print('window: ', window, window + 5, len(layer))
                    #layer_array = np.sort( layer[window:window+5] ).tolist()
                    #layer_win = tuple( layer_array )
                    """
                    if random_meas == {}:
                        random_meas[layer_win] = 1 ## This is only for the first iteration...
                    else:
                        rm_lays = list( random_meas.keys() )
                        for i in rm_lays:
                            l_targ = layer_array
                            l_key  = i
                            
                            if l_key == l_targ:
                                random_meas[layer_win] += 1
                                break
                            else:
                                random_meas[layer_win] = 1
                    """
                    layer = np.sort( layer[window:window+5] )
                    layer = tuple( layer )
                    #print('check layers: ', layer[:5])
                    
                    if layer in random_meas:
                        random_meas[layer] += 1
                    else:
                        random_meas[layer] = 1


            haar = self.Haar_meas( random_meas )
            window_Hmeas.append( haar )
            window_Rmeas.append( random_meas )
                
        #print( random_meas )
        return window_Rmeas, window_Hmeas
 
    def Expressibility_Measure( self, randencoding, dev ):
        """
        This is going to have to be injected like a prob wihtin the backend at 
        different points of building the circuit
        
        Note the hard part here, just a little thinking, is the random sampling of the circuit
        
        Have to figure a way to interpret this metric, not too bad though
        """
        express = []
        rand_meas, haar_meas = self.random_meas( encoding_method=encoding )
        for index, window in enumerate( rand_meas ):
            distr1 = {k:v/sum(rand_meas[index].values()) for (k,v) in rand_meas[index].items()}
            distr2 = {k:v/sum(haar_meas[index].values()) for (k,v) in haar_meas[index].items()}
            #print(distr1)
            expr = self.KLD( distr1, distr2 )
            express.append( expr )
 
        return express
    
    
    def KLD( self, random_meas, haar_meas ) -> float:
        """
        This will just be the formula for KLD of two distributions
        Note:The information gain when distr.2 (database) is used instead
             of distr.1 (input)
        """
        D_kl = []
        states = random_meas.keys()
            
        for i in states:
            D_kl.append(random_meas[i]*np.log(random_meas[i]/haar_meas[i]))
            
        KLD = np.sum( D_kl )
        return KLD

def main_analysis():
    # Params: num_bands, timbre_qubits, spine_qubits, qubits, distr1: list, database_states: list
    num_bands = 3
    
    spine_qubits = 3
    timbre_qubits = 3
    qubits = timbre_qubits + spine_qubits
    
    # s = Expressibility( num_bands, timbre_qubits, spine_qubits, qubits, [], qstates )
    # express_s, pk_seqs = s.Expressibility_Measure( 'static', 'ibmq_manhattan' )
        
    p  = Expressibility( num_bands, timbre_qubits, spine_qubits, qubits, [], [] )
    distr_1, database_bundle = p.random_meas( M=50, encoding_method='pkbse', big_circ=True, dev='ibmq_montreal' )
    rand, haar = p.Express_larger( distr_1, 5, database_bundle )
    
    return [rand, haar]

## Note: also try to get analyze the population of samples and sparseness between layers!
##----------------------------------------------------------------------------
if __name__ == "__main__":
    ## Test the Haar_measure:
    qc = QuantumCircuit( 1 )
    
    expr = Expressibility( qc, 3, 3, 1, {}, {} )
    haar_m = expr.random_measure()
    
    plot_histogram(haar_m); plt.show()
    
    print('testing: ', haar_m)
    
    
    