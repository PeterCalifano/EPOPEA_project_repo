# EPOPEA_project_repo
Code repository of the Enceladus Plume and Ocean Prospecting and Exo-Astrobiology (EPOPEA) Phase A analysis: 

To clone this repository:
1) Install git for Windows (64bit version)
2) Find the directory where you want the cloned repository to be
3) Right click on any point --> Git BASH here --> Make sure the prompt is opened in the current directory
4) Write "git clone https://github.com/PeterCalifano/EPOPEA_project_repo.git" and wait. Download of the repository should start
5) At the end, check that a new folder named as the repository has been created and that it contains all the files

In MATLAB work in the folder of the repository to use git integration. You will immediately see the difference with respect to an usual folder.
Legend of the symbols in the "git" column for a file: 
1) Green circle: Remote and Local origin have the same version 
2) Blue square: Local origin has one or more modification (even spaces count)
3) White circle: Local file only. Remote origin does not contain it: not synched.
4) Red circle with ! inside: conflict between Remote and Local file during push or pull procedure. Automatic merge has failed.

Basic commands from MATLAB GUI (and Git Bash):
1) Right click (RC) Source control (SC) --> Add to git to add a file with the white circle (git add <file>)
2) RC --> SC --> View and commit changes (+ write MEANINGFUL comment. Use capital letter as first letter because aesthetics matter. Lol) to issue a change from Local to      Remote origin (git commit)
3) RC --> SC --> Push, to forward the update to the Remote origin synching it with the Local stage status (git push)
4) RC --> SC --> Pull, to synch Local origin if Remote origin is ahead (git pull)
  
For Git Bash only: git status to check the Local origin status
Git GUI is also available by RC and can do everything descrived above but in a way similar to git integration for MATLAB.


