# Automated Crystal Orientation Mapping in py4DSTEM using Sparse Correlation Matching

by Colin Ophus, Steven E Zeltmann, Alexandra Bruefach, Alexander Rakowski, Benjamin H Savitzky, Andrew M Minor, and MC Scott

## Abstract

Crystalline materials used in technological applications are often complex assemblies composed of multiple phases and differently oriented grains. Robust identification of the phases and orientation relationships from these samples is crucial, but the information extracted from the diffraction condition probed by an electron beam is often incomplete. We therefore have developed an automated crystal orientation mapping (ACOM) procedure which uses a converged electron probe to collect diffraction patterns from multiple locations across a complex sample. We provide an algorithm to determine the orientation of each diffraction pattern based on a fast sparse correlation method. We test the speed and accuracy of our method by indexing diffraction patterns generated using both kinematical and dynamical simulations. We have also measured orientation maps from an experimental dataset consisting of a complex polycrystalline twisted helical AuAgPd nanowire. From these maps we identify twin planes between adjacent grains, which may be responsible for the twisted helical structure. All of our methods are made freely available as open source code, including tutorials which can be adapted to perform ACOM measurements on diffraction pattern datasets.

## Building the Paper

The paper has been translated from LaTeX to MyST Markdown, and can be rendered using the [Curvenote CLI (see install instructions)](https://curvenote.com/docs/cli/installing) and run:

```bash
curvenote start
```

The paper uses a [GitHub Action](https://curvenote.com/docs/web/github-action) to deploy the content to <https://rowanc1-py4dstem_paper.curve.space>

## Publishing Demonstration

Authors can write in Jupyter Notebooks using MyST markdown, or in the Curvenote platform, and then deploy to the Curvenote publishing infrastructure. The journal can aggregate content from multiple authors into a single user interface for browsing. Curvenote follows best practices for performance, accessibility, and archiving.
