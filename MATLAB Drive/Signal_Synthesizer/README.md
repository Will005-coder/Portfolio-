William Dakare: Lead Algorithm Developer
Noah Nir:  Data Manager
Samuel Yeboah-Asi: Data Visualization Specialist

1: File Dependencies:
 MATLAB Version: R2022b or later
 Toolbox Used: Signal Processing Toolbox; Statistics and Machine Learning Toolbox

2   a. To get results and conduct setup, open run `signalProcessingAlgorithm.m` in `./scripts.`. Make sure `helperScript.m` is in the same directory.
    b. To run tests on structure-object building `shouldTest()` and validating `helperScript.processSignals()`; pass parameters true in both 

3   We decided to create our directory in matlab drive. by doing this we
    were able to seamlessly collaborate with eachother. We structured the
    directory of our synthesizer by having all our code in the
    Signal_Synthesizer folder. This folder holds all our sub directories.
    Firstly in our main folder is our data folder. This Folder holds three
    subsections of clean, low noise, and noisy csv documents that provide the
    data for our scripts to use. Next we have the results folder. This folder
    holds all of our data figures and is clearly labeled in order to allow a
    outside reader to recognize and use our data. Finally is our scripts folder which holds all the scripts that we used in order to create and tranform our data into something that is usable.

4   Additionally, deliverables like `signalAnalysisResults.mat` and `signalAnalysisSummary.csv` are in the path `./results` with `./results/figure`

5   We decided to put our data in CSV format in order for it to be easily
    processed. We chose to have our first column as time in order to make the
    raw data legible and easy to use. We decided to name our signal columns by signal name followed by frequency as we found it to be the easiest to read when tyring to use the data of a particular HZ level and comparing it to different levels of noise. By doing this we found the data to be easy to read and more importantly seamless for the use of our code.



KEY FINDINGS:
- Noise misrepresents signal points, making it seem stronger at random intervals, thus reducing its quality.

- From previous analysis, the square wave is less robust than the sine and triangle waves, becasue it is more affected by noise.

- The moving average filter is very effective at reducing noise.
    With signal quality improvement rates of about 8% for robust wave-types like sine-wave
    And up to 16% for less robust waves like the triangle wave.

- Yes. Higher frequencies show more degradation than lower frequencies
    This is shown by 10 Hz having the highest degradation rate compared to 1 and 5 Hz.
    In another sense, the spontaneity of signal highs and lows in higher
    frequency signals makes their measured signals more prone to higher data misinterpretation



REFLECTION QUESTIONS
(Noah Nir): As Data Manager I found it incredibly important to create a
directory and data that was not only easy to read but usable in our code.By
doing this I found that the data was very easy to use by my partners and I
felt as if I was the base of our code.

(Samuel Yeboah-Asi: Visualization Specialist): Visualization is a huge part of projects like this; it helps draw out trends and is a catalyst in informed decision-making using data. For my part I made sure the core structure of all visualizations were balanced. All colors linewidth and naming were conventionally in sync. This ensures there's no headache for whoever utilizes the graphics made. Following the guide for visualizaion specialists also made selecting the appropriate graph type even easier.  

(William Dakare: Algorithm Developer): The core data organization method I used in my algorithm was the structure object. It was heavily commented on in the early parts to set the idea of its implementation, and benefits. I found using MATLAB's copilot effective in showing me relevant documentation. Beyond that using the "show/open" variable feature to get the content of a variable or structure was very useful in terms of back-tracing the points in my code that could have caused certain errors. My decision to extract the csv's into a structure object was for debugging benefits.

GROUP QUESTIONS
We found that working in a modular context was extremely important. This is because for one of us to do all the work it would take a much longer amount of time than dividing tasks into modules. This was especially helpful in algorithm development and data structures, because the smaller modules transfer low-risk of errors in the main code sequence. Additionally, we found that we learned from our past mistakes in creating code, as by properly labeling we were able to very easily read each other’s code and maintain a seamless production line. This closely mirrored real engineering work, as by working together efficiently we were able to make a cohesive project. It was also very important to keep the core algorithm simple for others to understand.
