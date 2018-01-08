**mindthegap**
=============

mindthegap is an opensource tracing software developed as a part of a project in **Google Summer of Code 2014**. It vectorizes bitmaps (currenlty supports png only) without introducing gaps or overlaps between adjacent areas.

- **Authors**: *Dhruv Kohli*, *Dr. Rembrandt Bakker*, *Dr. Piotr Majka*
- **Mentoring Organization**: *International Neuroinformatics Coordinating Facility*
- [Link](https://github.com/chiggum/Vectorization-of-brain-atlases) to **GSoC'14 project**: *Real-time Vectorization of brain atlases*
- [Link](https://scalablebrainatlas.incf.org/macaque/DB09) to a **Scalable Brain Atlas Template** where this tracing software is used.

Comparison
------------

**To the bottom of the page**

Usage
-------
- **For Linux**: `make -f makefile`. The binary will be generated in `bin` directory.
- **For windows**: Directly use exe file in `exe` directory.

Use `mindthegap -h`, for further help.

**Note**
-----
If the input image is such that different regions of the image can have multiple shades of same color then use the noisy switch by `-z`.

Components of overall algorithm
---------------------------
- Popping out boundaries between different colored regions (refer [documentation](https://chiggum.github.io/mindthegap/docs/documentation.pdf))
- Search Algorithm (DFS)
- Connected Component Labelling
- Dangerous connections removal (refer [documentation](https://chiggum.github.io/mindthegap/docs/documentation.pdf))
- Dissolving regions (refer [documentation](https://chiggum.github.io/mindthegap/docs/documentation.pdf))
- Bezier curve fitting over digitized curves **[1]**
- Posterization and Median Blurring (noisy switch)

References
-----------
1. Philip J. Schneider. “**An Algorithm for Automatically Fitting Digitized Curves**”. In Graphics Gems, Academic Press, 1990, pp. 612–626.


Comparsion
-----------
A comparison between currenlty most popular open source tracing software, AutoTrace and our tracing software, mindthegap.


**Input Bitmap**:

<img src="https://chiggum.github.io/mindthegap/docs/atlas_219.png" alt="Input bitmap image" width="300" height="600"/>

**AutoTrace SVG**:

<img src="https://chiggum.github.io/mindthegap/docs/output.svg" alt="Autotrace output" width="400" height="600"/>

**mindthegap SVG**:

<img src="https://chiggum.github.io/mindthegap/docs/mindthegap.svg" alt="mindthegap output" width="400" height="600"/>