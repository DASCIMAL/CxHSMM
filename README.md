# CxHSMM - The Coxian Hidden Semi-Markov Models
### Source codes for the Coxian Hidden Semi-Markov Models (CxHSMM) [Matlab/Octave]

***Copyright by Dinh Phung, Thi Duong and Hung Bui, 2005***
Version: 1.0

This package implements the **Coxian Hidden Semi-Markov Models** (CxHSMM) described in the following papers:

Duong, T., Phung, D., Bui, H. and Venkatesh, S, *Efficient Coxian Duration Modelling for Activity Recognition in Smart Environment with the Hidden semi-Markov Model*,  In Second International Conference on Intelligent Sensors, Sensor Networks and Information Processing, Melbourne, 5-6 December 2005.

Duong, T., Phung, D., Bui, H. and Venkatesh, *Efficient duration and hierarchical modeling for human activity recognition,* Artificial Intelligence (AIJ), 173(7-8):830-856, 2009.

##### Other related papers (e.g., use this code for baseline comparison, etc) include:

Phung, D., Duong, T., Bui, H. and Venkatesh, *Topic Transition Detection Using Hierarchical Hidden Markov and Semi-Markov Models,  In ACM Int. Conf on Multimedia (ACM-MM), Singapore, 6--11 Nov. 2005.

Duong, T., Bui, H., Phung, D. and Venkatesh, *Activity Recognition and Abnormality Detection with the Switching Hidden Semi-Markov Model*,  In IEEE Int. Conf. on Computer Vision and Pattern Recognition (CVPR), pages 838-845, San Diego, 20-26 June 2005.


**Disclaimer**: We have made our best effort in ensuring fairness in acknowledging existing codes and any materials we used. However, if you have any question/concern, please write to us.

### Using the package

**run (and look at) the following scripts**
```
- demo1_setup
- demo1_learn
- demo1_classify
```
  
  
### Citation information
  If you use the codes for your work, please use the following citation information.


```@INPROCEEDINGS { duong_phung_bui_venkatesh_issnips05,
    TITLE = { Efficient Coxian Duration Modelling for Activity Recognition in Smart Environment with the Hidden semi-Markov Model },
    AUTHOR = { Duong, T. and Phung, D. and Bui, H. and Venkatesh, S. },
    BOOKTITLE = { Second International Conference on Intelligent Sensors, Sensor Networks and Information Processing },
    YEAR = { 2005 },
    ADDRESS = { Melbourne },
    MONTH = { 5-6 December },
    ABSTRACT = { In this paper, we exploit the discrete Coxian distribution and propose a novel form of stochastic model, termed as the Coxian hidden semi-Makov model (Cox-HSMM), and apply it to the task of recognising activities of daily living (ADLs) in a smart house environment. The use of the Coxian has several advantages over traditional parameterization (e.g. multinomial or continuous distributions) including the low number of free parameters needed, its computational efficiency, and the existing of closed-form solution. To further enrich the model in real-world applications, we also address the problem of handling missing observation for the proposed Cox-HSMM. In the domain of ADLs, we emphasize the importance of the duration information and model it via the Cox-HSMM. Our experimental results have shown the superiority of the Cox-HSMM in all cases when compared with the standard HMM. Our results have further shown that outstanding recognition accuracy can be achieved with relatively low number of phases required in the Coxian, thus making the Cox-HSMM particularly suitable in recognizing ADLs whose movement trajectories are typically very long in nature. },
     TIMESTAMP = { 2010.08.11 },
}
```
```
@ARTICLE { duong_phung_bui_venkatesh_aij09,
    TITLE = { Efficient duration and hierarchical modeling for human activity recognition },
    AUTHOR = { Duong, T. and Phung, D. and Bui, H. and Venkatesh, S. },
    JOURNAL = { Artificial Intelligence (AIJ) },
    YEAR = { 2009 },
    NUMBER = { 7-8 },
    PAGES = { 830--856 },
    VOLUME = { 173 },
    ABSTRACT = { A challenge in building pervasive and smart spaces is to learn and recognize human activities of daily living (ADLs). In this paper, we address this problem and argue that in dealing with ADLs, it is beneficial to exploit both their typical duration patterns and inherent hierarchical structures. We exploit efficient duration modeling using the novel Coxian distribution to form the Coxian hidden semi-Markov model (CxHSMM) and apply it to the problem of learning and recognizing ADLs with complex temporal dependencies. The Coxian duration model has several advantages over existing duration parameterization using multinomial or exponential family distributions, including its denseness in the space of non-negative distributions, low number of parameters, computational efficiency and the existence of closed-form estimation solutions. Further we combine both hierarchical and duration extensions of the hidden Markov model (HMM) to form the novel switching hidden semi-Markov model (SHSMM), and empirically compare its performance with existing models. The model can learn what an occupant normally does during the day from unsegmented training data and then perform online activity classification, segmentation and abnormality detection. Experimental results show that Coxian modeling outperform a range of baseline models for the task of activity segmentation. We also achieve a recognition accuracy competitive to the current state-of-the-art multinomial duration model, whilst gain a significant reduction in computation. Furthermore, cross-validation model selection on the number of phases K in the Coxian indicates that only a small K is required to achieve the optimal performance. Finally, our models are further tested in a more challenging setting in which the tracking is often lost and the set of activities considerably overlap. With a small amount of labels supplied during training in a partially supervised learning mode, our models are again able to deliver reliable performance, again with a small number of phases, making our proposed framework an attractive choice for activity modeling. },
    PUBLISHER = { Elsevier },
    TIMESTAMP = { 2010.08.11 },
}
```
