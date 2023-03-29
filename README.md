# Hemolymph viscosity in hawkmoths and its implications for hovering flight - Dataset:
---

Dataset contains the average values, with their standard deviations, for each species of hawkmoths we probed in the experiment. All variables are contained in the same file. These data were used to run all analyses contained in the manuscript.

Hemolymph viscosity was obtained by tracking nickel nanorods in droplets of hemolymph. All calculations and how we did it is contained in the main paper. The intermuscle distances, on the other hand, were calculated in micro-CT images obtained for each species. Given that we are using the micro-CT scans in other research questions, we are not uploading them together with the dataset. However, they are available upon request.  

We are also uploading the R codes and the phylogenetic tree used to run all analyses.
Phylogenetic tree was taken from: 
Kawahara, A. Y., & Barber, J. R. (2015). Tempo and mode of antibat ultrasound production and sonar jamming in the diverse hawkmoth radiation. Proceedings of the National Academy of Sciences, 112(20), 6407-6412.

## Description of the Data and file structure

We are providing two files: the dataset and the codes. 
The metadata of the data is as follows:

METADATA OF hemolymph_data.csv

In the columns we have the variables, in rows we have the individuals. We obtained mean values by averaging the individuals of the species.

COLUMN A: species - the hawkmoth species <br>
COLUMN B: m.visc - mean viscosity for the species, in mPa*s. We measured viscosity using magnetic rotational spectroscopy (MRS). For a full description of how it was measured, please check the main manuscript. <br>
COLUMN C: visc.std - standard deviation of the viscosity for the species <br>
COLUMN D: body.length - mean distance between the distal end of the abdomen and the distal end of the head. Unit: cm <br>
COLUMN E: length.sd - standard deviation of the distance between the distal end of the abdomen and the distal end of the head. <br>
COLUMN F: body.width - mean distance between the base of the wings. Unit: cm. <br>
COLUMN G: width.sd - standard deviation of the distance between the base of the wings. <br>
COLUMN H: volume - mean inferred body volume. We inferred body volume using the equation of a circular cylinder using body length and width measurements. Unit: cm cubed <br>
COLUMN I: vol.sd - standard deviation of the inferred body volume. <br>
COLUMN J: prob.length - mean length of the proboscis. Measured from the head to the tip. Unit: cm. <br>
COLUMN K: forewing.length - mean distance between the base of the wing to the outermost tip. Unit: cm. <br>
COLUMN L: forewing.sd - standard deviation of the distance between the base of the wing to the outermost tip. <br>
COLUMN M: forewing.max - highest measurement within a species of the distance between the base of the wing to the outermost tip. Unit: cm <br>
COLUMN N: N - number of individuals probed for hemolymph viscosity. <br>
COLUMN O: sp.in.tree - name of the species in the phylogenetic tree. <br>
COLUMN P: group - subfamily of Sphingidae to which the species belongs. <br>
COLUMN Q: label - just a numbered list to add labels in the figures. <br>
COLUMN R: mean.slit - average distance between the halves of the longitudinal flight muscles. Used to infer the pore structure of the thorax. The distances were measured in five equidistant places along the thorax of the hawkmoths. Unit: mm. <br>
COLUMN S: sd-slit - standard deviation of the distance between the halves of the longitudinal flight muscles. <br>
COLUMN T: maximum - maximum distance between the halves of the longitudinal flight muscles. Used to infer the pore structure of the thorax. The distances were measured in five equidistant places along the thorax of the hawkmoths. Unit: mm. <br>
COLUMN U: minimum - minimum distance between the halves of the longitudinal flight muscles. Used to infer the pore structure of the thorax. The distances were measured in five equidistant places along the thorax of the hawkmoths. Unit: mm. <br>

The code was run in R software v4.2.2.
Packages used:
ape - v5.6.2
phylolm - v2.6.2
phytools - v1.2.0
geiger - v2.0.10
treeio - v1.22.0
RColorBrewer - v1.1.3
scales - v1.2.1

Phylogenetic tree was taken from:
Kawahara, A. Y., & Barber, J. R. (2015). Tempo and mode of antibat ultrasound production and sonar jamming in the diverse hawkmoth radiation. Proceedings of the National Academy of Sciences, 112(20), 6407-6412.


## Sharing/access Information

The file structure and files can be seen and downloaded from:
https://github.com/alexandrepalaoro/hemo-viscosity 

The phylogenetic tree was taken from another paper:
Kawahara, A. Y., & Barber, J. R. (2015). Tempo and mode of antibat ultrasound production and sonar jamming in the diverse hawkmoth radiation. Proceedings of the National Academy of Sciences, 112(20), 6407-6412.
