Instructions for running locally
1.Move to build directory
2.cmake ..
3.make 








Instrucitons for running on Augsburg cluster
1.# To login to your account in Augsburg (Remember the passphrase in your home folder)
ssh -i rahul_kumar_key.pem rahul@cluster.math.uni-augsburg.de


2. # passphrase
SEmGU8v4

3. # To copy from Desktop to Augsburg cluster

scp -r /home/rkumar/Desktop/Collaborations/UH/tvf_c++_new/ErrorIPDG.ufl  rahul@cluster.math.uni-augsburg.de:/homes/gast/rahul


4. # To copy from Augsburg cluster to Desktop 

scp -r rahul@cluster.math.uni-augsburg.de:/homes/gast/rahul/ErrorIPDG1.ufl /home/rkumar/Desktop/Collaborations/UH/tvf_c++_new 

