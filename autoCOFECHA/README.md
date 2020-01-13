----
Automate COFECHA.exe README
----

The COFECHA program was not created by me but can be found at <https://www.ltrr.arizona.edu/pub/dpl/>. [The Laboratory of Tree Ring Research at University of Arizona](https://ltrr.arizona.edu/resources) provides a list of other dendrochronology software.

This is a bash script that you can use to automate COFECHA. The bash script contains two nested for loops: one for the files you want to pull up and the other to run through multiple splines. The GoBaby script calls COFECHA and inputs the necessary answers. (Shout out to this [Stack Overflow answer](https://stackoverflow.com/questions/41883729/how-to-run-exe-with-inputs-on-terminal-in-unix-linux) for the template). If you need to change something when talking to COFECHA, you would edit GoBaby. 
