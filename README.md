# Registration
This repository registers two images using maximization of mutual information.\
The `Mutual_Information_Computation` function computes mutual information using rigid body registration. It takes as input the angle of rotation and row-wise and column-wise translation.
`Mutual_Information_Computation` is a vectorized function. It does not contain any loop. Vectorization is used in order to increase speed of computation.\
The formula for mutual information contains entropy and joint entropy. Entropy and joint entropy are computed using the standard formulae, as well as our proposed formulae implemented in `Efficient-Entropy`.
In the function, both the formulae are implemented. They are demarcated and commented separately. While running, any one set of formulae should be uncommented, while the other set should be commented out.\
This work has not been published anywhere. It has been implemented in SimpleITK using PET-CT and MR images of the brain, and is under revision for possible publication.
