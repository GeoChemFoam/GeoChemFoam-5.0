These 8 files are the 8 batch scripts that constitute the steps of the workflow which 
will automatically run tutorial 1 in full, where each job submits the next, in order. 
The scripts use gcfoam module.

To launch this workflow, 

1 change the username and project code to your username and project code, 
2 copy all 8 files into the parent (Ketton) directory above
3 change directory to the parent direcorty: cd ..
4 launch the workflow but submitting the first batch script: sbatch workflow_step_1.bat

