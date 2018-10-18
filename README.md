# Acoustic-Beamforming
Various beamforming files written in MATLAB. Allows for convection, different types of steering vectors and more. They are known to work with R2016b (9.1.0.441655).

Requires time pressure signal of microphones as a matrix representation.
Use developCSM to construct CSM with a given sampling frequency.
Use any beamform algorithm to perform beamforming.

Q: "I just want to beamform...?"

A: If you have your time pressure signals and the microphone configuration you should be fine with only using <developCSM.m> and <FastBeamforming3.m>. An example can be seen in <Main1_Basic.m> in my 'Beamforming-Simulations' repository.

Q: "But I don't have array data...:("

A: You can simulate array data using <simulateArraydata.m>. It requires the source position, microphone positions, sample frequency and other parameters. See the function description. Allows for multiple sources and flow. The output is time pressure signals for the different microphones. This is then used as an input to <developCSM.m>. An example where it is used can be seen in <Main1_Basic.m> in my 'Beamforming-Simulations' repository.

Q: "What is going on in <FastBeamformingX.m>?"

A: The beamforming algorithm is rewritten to apply fast vector manipulation. This is harder to comprehend (we go from matrices to vectors and do beamforming). To understand what is going on, you can use <SarradjBeamforming.m> which is exactly the same given the steering vector, e.g. 'form3' corresponds to <FastBeamforming3.m>. <SarradjBeamforming.m> only considers a single scan point, i.e. you have to build your own loop outside the function to evaluate a scan plane. In the end it is easier to understand 'what is going on', but it will also be severely slower than the other function.

Q: "What is <..Conv.m>?

A: It just includes another parameter for the flow speed, i.e. it corresponds to the same script without 'Conv' at the end of the name for U = [0 0 0]. Files will be removed/re-arranged later. (This scripts were used as seperate to investigate performance between them, i.e. some are now obsolete as you can simply put a zero vector in 'Conv'.)

Q: "What is this <..mod.m>?"

A: It is a small modification to have the steering vectors exactly as how Sijtsma describes it in his tech papers. That is, the powers of your scan points will always be relative to a reference distance. See <FastBeamforming3mod.m> comments for additional info. Sijtsma's steering vector can also be called from <SarradjBeamforming.m> and using 'pieter' as steering option.

Q: "I would love to experiment with this beamforming stuff... if only I could use my mouse to click for multiple source positions on a map, do the same to construct my microphone array, and see the result with beamforming...!"

A: Are you lucky or what? By chance (ahem) I made these kind of scripts! It can be found in... (continue answer and add the files).
