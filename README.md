In this project, we, the students, were tasked with specializing in a selected area of research based on scientific method, e.g., collecting supplementary information based on literature search and practical and other sources and combining it with own knowledge into a project report. Furthermore, we had to complete a comprehensive independent project, including drawing up a project plan with milestones, reporting partial results and progress, and writing a standard project report. The project was a part of the 5th year Cybernetics and Robotics curriculum.

# Abstract

Sepsis is one of the leading causes of death for hospitalized patients. Early symptoms are challenging to identify, as they are subtle, and other conditions with a common physiologic response meet the same criteria of diagnosis. During sepsis, the depressed vascular properties of the circulatory system are the main contributor to decreased blood perfusion. The reduced vascular tone will cripple the autoregulatory mechanisms that the body controls to maintain appropriate blood perfusion by controlling pressure and vascular dilation. 


In a healthy patient, the autoregulation is observed as distinct slow variations in the hemodynamic parameters, such as systemic vascular resistance (SVR) and arterial compliance. Provided velocity and pressure measurements, we can calculate the arterial impedance, which can be analogized to a 2-element Windkessel model (WK). The model resembles a parallel RC circuit, in which the resistance is the SVR, and the capacitance is the arterial compliance. During septic shock, the autoregulatory mechanisms will be unstable and show a more stochastic behavior in the 2-element WK model as it attempts to stabilize the blood flow. This stochastic behavior should be derivable as an increase in the magnitude of the low autoregulatory frequencies and its neighbors. 


The 2-element WK model parameters were found from 3-minute blood pressure and velocity measurements of seven patients. The parameters were then transformed into the frequency domain through a digital Fourier transform (DFT), where the possibility of distinguishing the severity of the septic shock from looking at the average magnitude of the working frequencies of the autoregulation was investigated. Although there seems to be some correlation, there is not enough evidence to employ this as a diagnostic measurement for sepsis. I do believe there is potential for further and more comprehensive studies on this matter, and I do see the potential to develop related research with this project's results and potential for development as a baseline.

# Repocitory Outline

The **report** folder holds the pdf and latex source code for the project report, while the **scripts** folder contains all the scripts used while working on the project.

The scripts are kind of messy, as they are immediate and straightforward scripts to analyze the datasets and not a production-ready application. I have not bothered to clean up the code too much or add any excessive commentary, as it is uploaded for source control during the project.


