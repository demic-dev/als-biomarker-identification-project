import xml.etree.ElementTree as et
import matplotlib.pyplot as plt # type: ignore
import pandas as pd # type: ignore
import numpy as np
import re

from modules import Annotations, Samples

class Sequencing:
    #region Default Methods

    def __init__(self, samples: Samples.Samples, annotations = Annotations.Annotations) -> None:
        self.__samples = samples
        self.__annotations = annotations

        return

    #endregion

    #region Setters and Getters
    
    def get_samples(self) -> pd.DataFrame:
        return self.__samples
    
    def get_annotations(self) -> pd.DataFrame:
        return self.__annotations
    
    #endregion
  
    #region Statistical Analysis

    #endregion

    #region Sample Description
    
    def get_sample_distribution_by_column(self, column: str) -> plt.Axes:
        return self.__sample_annotations[column].value_counts().plot(kind='bar', rot=45)
    
    #endregion
