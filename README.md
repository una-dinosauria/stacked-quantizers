Stacked Quantizers
==================

This is the code for the paper

Julieta Martinez, Holger H. Hoos and James J. Little: *Stacked Quantizers for Compositional Vector Compression*, available at: http://arxiv.org/abs/1411.2173

This code was mostly written by [Julieta Martinez](https://github.com/una-dinosauria/).

### Datasets

The demo requires you to download the `convnet1m-128` dataset, and put it into a `/data` folder at the top directory.
* The datasets for training, queries and database are available in:  [features_m_128.mat](https://drive.google.com/file/d/0BxWzojlLp259aUt4Z1ZCMzVDQlk/view?usp=sharing)
* The ground truth is available in: [features_m_128_gt.mat](https://drive.google.com/file/d/0BxWzojlLp259cmJkeG40S3oxR28/view?usp=sharing)

`SIFT1M` and `GIST1M` can be downloaded from [INRIA](http://corpus-texmex.irisa.fr/).

### Demo

For a demonstration of approximate nearest neighbour search 

1. Compile search utilities by running `compile.m` from the top directory.
2. Run `demo.m` from the top directory as well.

### Citation

If you find this code useful, please consider citing our paper:

Julieta Martinez, Holger H. Hoos, and James J. Little.: *Stacked Quantizers for Compositional Vector Compression.* arXiv preprint arXiv:1411.2173 (2014).

### Acknowledgements

* The code under `/OPQ` was taken from [Cartesian k-means](https://github.com/norouzi/ckmeans) by [Mohammad Norouzi](https://github.com/norouzi)
* The code that draws the beautiful `recall@N` plots under `/SQ/util/eval_recall_vs_sel.m` is by [Kaiming He](http://research.microsoft.com/en-us/um/people/kahe/)



