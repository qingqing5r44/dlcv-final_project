DL for CV final project:      

3D human model reconstruction from a single RGB image         

The 3D human body reconstruction process mainly consists of three parts (OpenPose, SMPLify-X, PanoMan) and five actions. They are two-dimensional key point extraction, estimation of SMPL model parameters, acquisition of three-dimensional joint points, PanoMan model parameter initialization and reconstruction optimization. First, use the input RGB image to extract 2D key points, estimate SMPL parameters and acquire reference values of 3D joint data. Then use the estimated values of SMPL parameters to initialize the PanoMan model shape and motion mixed base. Finally, the joint data reference value is used to constrain the joint points of PanoMan to optimize the model results, so as to obtain a model closer to the human body in the picture.                 

OpenPose folder includes the code for reimplementing the OpenPose network, training and testing the trained model.                 

SMPLify-X folder includes the official code provided by the team. Here I made small changes to fit my need.            

PanoMan folder includes the code of implementing the PanoMan human model, the trained MLP for transferring the SMPL pose parameters into PanoMan mixed coefficients, reconstruction procedure, and the qt ui file and scripts for the 3D reconstruction system.                 
