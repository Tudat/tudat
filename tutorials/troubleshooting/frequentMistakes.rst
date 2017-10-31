.. _troubleshootingFrequentMistakes:

Frequently made  coding mistakes
================================
This section sums op some of the most made mistakes leading to failures either in runtime or during compilation. This page is under construction and will be elaborated in the future.

- Accessing vector elements out of vector bounds::
      
    Eigen::Vector2d vectorOfLengthTwo(2.0, 3.0);
    double intermediateResult = vectorOfLengthTwo(2); 

 This will actually work, but the value obtained ``intermediateResult`` is the value located just next to the ``vectorOfLengthTwo`` in your memory and could be anything. 

