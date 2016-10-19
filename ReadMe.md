**mindthegap**
=============

mindthegap is an opensource tracing software developed as a part of a project in **Google Summer of Code 2014**. It vectorizes bitmaps (currenlty supports png only) without introducing gaps or overlaps between adjacent areas.

- **Authors**: *Dhruv Kohli*, *Dr. Rembrandt Bakker*, *Dr. Piotr Majka*
- **Mentoring Organization**: *International Neuroinformatics Coordinating Facility*
- [Link](https://github.com/chiggum/Vectorization-of-brain-atlases) to **GSoC'14 project**: *Real-time Vectorization of brain atlases*
- [Link](https://scalablebrainatlas.incf.org/macaque/DB09) to a **Scalable Brain Atlas Template** where this tracing software is used.

Comparison
------------
A comparison between currenlty most popular open source tracing software, AutoTrace and our tracing software, mindthegap. **Right click and open in new tab/window to see the difference in outputs of two tracing softwares**.

<img src="https://chiggum.github.io/mindthegap/docs/atlas_219.png" alt="Input bitmap image" style="width: 100px;"/>
<img src="https://chiggum.github.io/mindthegap/docs/output.svg" alt="Autotrace output" style="width: 100px;"/>
<img src="https://chiggum.github.io/mindthegap/docs/mindthegap.svg" alt="mindthegap output" style="width: 100px;"/>

Usage
-------
- For Linux: `make -f makefile`. The binary will be generated in `bin` directory.
- For windows: Directly use exe file in `exe` directory.

Use `mindthegap -h`, for further help.

Note
-----
If the input image is such that different regions of the image can have multiple shades of same color then use the noisy switch by `-z`.