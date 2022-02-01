# QuiKo-v2.0
This is a quantum beat generator that constructs a set of beats and rhythmic patterns based on the audio features extracted from an input audio file. The beat is constructed from a a database of audio and drum samples that reflect the features of the 
QuiKo Workflow and Walk Through
SET UP
Input Signal: 1st measure of HER - Feel a way https://drive.google.com/file/d/1ROxWiFnVbZDsTWD3Ztx1_dCWfv4ucPbt/view?usp=sharing
 as .wav, mp3 etc. I have the audio file. The first measure should be the first two seconds of the song.

Full Song: (https://www.youtube.com/watch?v=hvqQo87bdys) 
Signal Prep:

Apply low pass filter to input signal call the resulting signal ‘her_low.wav’
Apply band pass filter to input signal call the resulting signal ‘her_mid.wav’
Apply high pass filter to input signal call the resulting signal ‘her_high.wav’

      We now have three filtered versions of the same input signal. This will be used in the 

Drum Sample Database Preparation:

Gather a collection of drum samples and audio samples to be used in the generated beat. (just put them in a folder labeled ‘Audio Files’)

Apply a low, band pass and high pass filter to each sample in the collection as we did with the input signal. Keep track of the filtered versions of the drum samples. There should be a low, mid and high version of each sample.

For each of the sample’s filtered versions Apply the HPSS algorithm to get the percussive and harmonic parts. 

For the percussive part find the values of the peaks (spikes) in the signal, sum their values together then divide by the max peak (the value of the maximum peak) → This will be our  values for the U3-gate



For the harmonic part take the FFT and find the highest 3 peaks in the spectrum and take their weighted average. → This will be our 

Also for the harmonic part calculate it spectral centroid → This will be our  value.

The calculations in (3) will be done for each of the filtered versions of the original samples. The values will be mapped to specific qubits. The angles calculated in (3) for the low version of the sample will be mapped to q0 using a U3 gate in qiskit. The angles for the mid  will be mapped to q1, and high mapped to q2

		Low angles  →  q0
		Mid angle     →  q1 
      		High angle   →  q2

Note: Also note that the qubits are put in superposition first before we apply the U3 gates. 


Measure this circuit on a simulator or on a real device. I currently have only run this section on the simulator. This circuit should be executed, with 2048 shots, for each unique sample in the database with their specific angles obtained from the feature extraction in (3) and (4). 

From here there should be a probability distribution generated for each sample, and should be saved to a dictionary. This will be used later to determine which samples get mapped to a specific subdivision.

SIGNAL PREP

Just like in the previous section for preparing the database we will extract features from the input signal using the filtered versions of the input signal that we prepared in the beginning (signal prep). 

Using the original input signal (unfiltered) estimate the tempo using the librosa.beat.tempo function in librosa. This is given in BPM so calculate how many beat per second bps. With this you can determine how many seconds per measure if you also know the time signature. I am currently inputting the time signature manually.

Also determine how subdivisions per measure. This is based off how many qubits you wish to allocate for encoding subdivision information (a sort of qubit budget). In my case I allocated 3 qubits resulting in 8 subdivisions per measure.

For each of the filtered versions of the input signal apply the HPSS algorithm to get the percussive and harmonic parts. We now want to do feature extract as in the last section but here we want to do this for each of the 8 subdivisions for both the percussive and harmonic parts. So now divide the percussive and harmonic parts into 8 segments. Each segment will be associated with a binary code



Perform the same calculations from (3) and (4) from the previous section (Drum Sample Database Preparation) for each of the 8 segments.







QUANTUM CIRCUIT

After the feature extract is complete we need to assemble the quantum circuit below. We have three registers. 

b0, b1, b2 → This is the “Timbre” register which the U3gates with angle parameters calculated in last section will be applied. 

q0, q1, q1 → This register will represent which subdivision a specific sample occurs.
The last register (last three qubits in the circuit) will be cor measurement of the first register in (1).


Figure 4: Quantum Circuit for input signal

Apply hadamard gates to the first two registers.

Apply the multi-controlled U gates (qiskit U3 gates in this case) entangling the a specific subdivision segment to the features extracted for those subdivisions. 

For example, for subdivision 101 (or 5 in this case) the ,, and  found for the filtered versions of the input signal, for this subdivision, U3gates with the corresponding parameters will be applied to qubits b0, b1, b2. 
This will be repeated for every subdivision.

As a result of (a) phase kick back occurs and affects the second register qubits (which we prepared in the hadamard basis). Because of this the value of the second qubit register has encoded with the overall phase of the unitary gates. This will specify the new subdivision for which the determined sample will occur. To extract this data we perform phase estimation on the second qubit register (implementing the inverse QFT) and then measuring the second register.

The block (b) in the diagram is what we will call the react oracle. It specifies a rule in how to change the code of the timbre. You can also think of it as specifying a rule in how the system should react depending on the state of the timbre register after block (a) in the diagram. 

This is up to the user to specify this section. In the diagram above we are specifying if the qubit (b0) representing the lower band is one and the the qubit (b1) representing the mid band then flip the value of the qubit representing the higher band. Or if he mid band qubit is one flip the high band qubit.

We then in block c entangle the first qubit register with the third register to measure its results for that measure on the third register. We do this because ideally we do not want to collapse the state of the first qubit register as they will be the initial step of the next measure. 

INITIAL RESULTS

As a refresher the input signal was the 1st measure of HER - Feel a way https://drive.google.com/file/d/1ROxWiFnVbZDsTWD3Ztx1_dCWfv4ucPbt/view?usp=sharing as .wav, mp3 etc. I have the audio file. The first measure should be the first two seconds of the song.

The Drum Sample DataBase used samples from the free 9th Wonder sound pack and samples from Splice (https://splice.com/home).

MEASUREMENT

The circuit in the previous section, with all the parameters for the input signal, was executed for 2048 shots using the simulator. A probability distribution was generated showing the probabilities of six bit binary strings. The first three bits in the binary string specifies the subdivision.

For each subdivision multiple timbre codes can occur (i.e. 000 | 001, 000 | 010, etc.). In addition, since we only have 8 states to assign, and if we had a database larger than 8 samples for more variety we would begin to have conflicting timbre codes. As a result, for each subdivision we compare the distribution of the timbre codes with the distribution in the sample database we previously prepared. This was done using Numpy correlate function to perform a cross correlation on the distributions. 



Figure 5: (Left) Shows the part of the binary string representing the subdivision (Right) Cross Correlation between the probability distribution of the quantum circuit for the input signal and the distribution for each of the samples in the database

After some experimentation the 5-10 min values for the correlate function produces the most rhythmic and varied Beat while the top 1-3 max values for the cross correlation produce a more overall atmospheric choice in sample. As a result multiple layers are produced for our new generated beat. Also note that the We could further organize and classify these outputs layers, but further investigation and research is required at this stage. 

RESULTS FROM SIMULATOR
Note: The outputs currently are a table showing the mapping of the samples to specific subdivisions.

TEST #1
Drum Sample Table Arrangement for the Max Cross_Correlation
 
Subdiv 1
Subdiv 2
Subdiv 3
Subdiv 4
Subdiv 5
Subdiv 6
Subdiv 7
Subdiv 8
Layer 1
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Layer 2
Samp_6.wav
Samp_6.wav
Samp_6.wav
Samp_6.wav
Samp_6.wav
Samp_6.wav
Samp_6.wav
Samp_6.wav


Drum Sample Table Arrangement for the Minimum Cross Correlation
 
Subdiv 1
Subdiv 2
Subdiv 3
Subdiv 4
Subdiv 5
Subdiv 6
Subdiv 7
Subdiv 8
Layer 1
Ban_K.wav
Samp_1.wav
Samp_1.wav
Ban_K.wav
Samp_14.wav
Samp_1.wav
Clk_H2.wav
Ban_K.wav
Layer 2
Samp_2.wav
Hi_Tom2.wav
Hi_Tom2.wav
Bub_K.wav
drum2.wav
Hi_Tom2.wav
Cd_H.wav
Samp_1.wav
Layer 3
Bub_K.wav
Samp_10.wav
Samp_10.wav
Conga_12.wav
Bar_K.wav
Samp_10.wav
Samp_15.wav
Samp_11.wav
Layer 4
Bar_K.wav
Samp_19.wav
Samp_19.wav
Crsh_K.wav
Samp_4.wav
Samp_19.wav
Samp_17.wav
Hi_Tom2.wav
Layer 5
Samp_4.wav
drum2.wav
drum2.wav
Samp_14.wav
Bngo_4.wav
drum2.wav
Bnt_Snr.wav
Conga_4.wav

Audio:  https://drive.google.com/file/d/12ua7QtLNR5jz5MtDrTitMlV7z-anlwLj/view?usp=sharing
TEST #2
Drum Sample Table Arrangement for the Minimum Cross Correlation (Did not find the layers for max correlation here)
 
Subdiv 1
Subdiv 2
Subdiv 3
Subdiv 4
Subdiv 5
Subdiv 6
Subdiv 7
Subdiv 8
Layer 1
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Layer 2
Bngo_4.wav
Bngo_4.wav
Bngo_4.wav
Bngo_4.wav
Bngo_4.wav
Bngo_4.wav
Bngo_4.wav
Bngo_4.wav
Layer 3
samp_4.wav
Bar_k.wav
Bar_k.wav
Bar_k.wav
Samp_19.wav
Samp_19.wav
Bar_K.wav
Samp_10.wav
Layer 4
Conga_8.wav
Conga_8.wav
Samp_4.wav
Samp_7.wav
Bar_k.wav
Conga_8.wav
Samp_10.wav
Cev_H1.wav
Layer 5
Samp_10.wav
Samp_4.wav
Samp_4.wav
Samp_16.wav
Samp_16.wav
Samp_4.wav
Samp_16.wav
Bmb_K.wav

Audio: https://drive.google.com/file/d/1LPRwOlp8QmoTW2P21TQTDJIMI2bkJHRL/view?usp=sharing


TEST #3
Drum Sample Table Arrangement for the Max Cross_Correlation
 
Subdiv 1
Subdiv 2
Subdiv 3
Subdiv 4
Subdiv 5
Subdiv 6
Subdiv 7
Subdiv 8
Layer 1
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Samp_16.wav
Layer 2
-
-
-
-
-
-
Samp_6.wav
Samp_6.wav


Drum Sample Table Arrangement for the Minimum Cross Correlation
 
Subdiv 1
Subdiv 2
Subdiv 3
Subdiv 4
Subdiv 5
Subdiv 6
Subdiv 7
Subdiv 8
Layer 1
Samp_2.wav
Samp_1.wav
Samp_1.wav
Conga_1.wav
Samp_4.wav
Samp_1.wav
Bmb_K.wav
Bngo_8.wav
Layer 2
Conga_1.wav
Hi_Tom2.wav
Hi_Tom2.wav
Ban_K.wav
Samp_2.wav
Hi_Tom2.wav
Samp_17.wav
Samp_14.wav
Layer 3
Samp_4.wav
Samp_10.wav
Samp_10.wav
Crsh_K.wav
Conga_4.wav
Samp_10.wav
Clk_H2.wav
Hi_Tom2.wav
Layer 4
Ban_K.wav
Samp_19.wav
Samp_19.wav
Bngo_8.wav
Samp_19.wav
Samp_19.wav
Dbn_K.wav
Samp_1.wav
Layer 5
Bub_K.wav
drum2.wav
drum2.wav
Bub_K.wav
Ban_K.wav
Samp_4.wav
Crs_Tmb_1.wav
Samp_13.wav

Audio:
https://drive.google.com/file/d/1g7cMS_GDzf1FVoOTnzvQk48nCh7V_553/view?usp=sharing

Note: Each of these are separate tests on the same input signal.





