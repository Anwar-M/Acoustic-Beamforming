# Acoustic-Beamforming
Various beamforming files written in MATLAB. Allows for convection, different types of steering vectors and more.

Requires time pressure signal of microphones as a matrix representation.
Use developCSM to construct CSM with a given sampling frequency.
Use any beamform algorithm to perform beamforming.

Q: "I just want to beamform...?"

A: If you have your time pressure signals and the microphone configuration you should be fine with only using <developCSM.m> and <FastBeamforming3.m>. An example can be seen in <Main1_Basic.m> in my 'Beamforming-Simulations' repository.

Q: "What is going on in <FastBeamformingX.m>?"

A: The beamforming algorithm is rewritten to apply fast vector manipulation. This is harder to comprehend (we go from matrices to vectors and do beamforming). To understand what is going on you can use <SarradjBeamforming.m> which is exactly the same given the steering vector, e.g. 'form3' corresponds to <FastBeamforming3.m>. <SarradjBeamforming.m> only considers a single scan point, i.e. you have to build your own loop outside the function to evaluate a scan plane. In the end it is easier to understand 'what is going on', but it will also be severely slower than the other function.

Q: "What is <..Conv.m>?

A: It just includes another parameter for the flow speed, i.e. it corresponds to the same script without 'Conv' for U = [0 0 0]. Files will be removed/re-arranged later. (This script was used as a seperate to investigate performance between different scripts, but is now obsolete.)

Q: "What is this <..mod.m>?"

A: It is a small modification to have the steering vectors exactly as how Sijtsma describes it in his tech papers. That is, the powers of your scan points will always be relative to a reference distance. See <FastBeamforming3mod.m> comments for additional info. Sijtsma's steering vector can also be called from <SarradjBeamforming.m> and using 'pieter' as steering option.
