# Acoustic-Beamforming
Various beamforming files written in MATLAB. Allows for convection, different types of steering vectors and more.

Requires time pressure signal of microphones as a matrix representation.
Use developCSM to construct CSM with a given sampling frequency.
Use any beamform algorithm to perform beamforming.

Q: "I just want to beamform..."
A: If you have your time pressure signals and the microphone configuration you should be fine with only using <developCSM.m> and <FastBeamforming3.m>. An example can be seen in <Main1_Basic.m> in my 'Beamforming-Simulations' repository.
