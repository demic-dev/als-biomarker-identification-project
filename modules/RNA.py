import xml.etree.ElementTree as et
import matplotlib.pyplot as plt # type: ignore
import pandas as pd # type: ignore
import numpy as np
import re
class RNASeq:
    """
    RNASeq class for processing RNA sequencing data.

    Attributes:
        verbose (bool): Verbosity mode. If True, prints out messages during data loading.
        length (int): Number of samples loaded.
        __genes_counts (DataFrame): DataFrame containing gene counts data.
        __sample_annotations (DataFrame): DataFrame containing sample annotations.

    Methods:
        parse_file(filename): Parses a data file and returns a DataFrame.
        load_single_data(filename): Loads gene counts data from a single file.
        load_bulk_data(filenames): Loads gene counts data from multiple files.
        load_annotations(file): Loads sample annotations from an XML file.
        get_std_genes(): Calculates standard deviation of gene counts across samples.
        get_median_genes(): Calculates median of gene counts across samples.
        get_mean_genes(): Calculates mean of gene counts across samples.
        get_genes_counts(): Returns the DataFrame containing gene counts data.
        get_sample_annotations(): Returns the DataFrame containing sample annotations.
    """

    #region Default Methods

    def __init__(self, verbose: bool = False) -> None:
        """
        Initializes RNASeq class.

        Args:
            verbose (bool, optional): If True, enables verbose mode. Defaults to False.
        """
        self.__genes_counts = pd.DataFrame()
        self.__sample_annotations = pd.DataFrame(columns = ["id","cns subregion", "tissue type", "sample group"])
        
        self.__index_label = "sequencing"

        self.length = 0

        self.__verbose = verbose
        return
    
    def __str__(self) -> str:
        return str(self.length)

    #endregion

    #region Setters and Getters
    
    def toggle_verbose(self) -> None:
        self.__verbose = not self.__verbose
        
    def get_verbose(self) -> bool:
        return self.__verbose
    
    def get_genes_counts(self) -> pd.DataFrame:
        """
        Returns the DataFrame containing gene counts data.

        Returns:
            DataFrame: Gene counts data.
        """
        return self.__genes_counts
    
    def get_sample_annotations(self) -> pd.DataFrame:
        """
        Returns the DataFrame containing sample annotations.

        Returns:
            DataFrame: Sample annotations.
        """
        return self.__sample_annotations
    
    #endregion
    
    #region Importing and utils
    
    def parse_file(self, filename: str) -> pd.DataFrame:
        """
        Parses a data file and returns a DataFrame.

        Args:
            filename (str): Path to the data file.

        Returns:
            DataFrame or None: DataFrame containing parsed data, or None if parsing failed.
        """
        try:
            df = pd.read_table(filename)
            sample_name = re.search("GSM\d+", filename).group()
            df.rename(index= df["gene/TE"], inplace=True)
            df.drop(columns=df.columns[0], axis=1, inplace=True)
            df.rename(columns={ df.columns[0]: sample_name }, inplace = True)
            return df
        except:
            return None
            
    def load_bulk_data(self, filenames: list[str], overwrite = True) -> None:
        """
        Loads gene counts data from multiple files.

        Args:
            filenames (list): List of paths to data files.
            overwrite (bool): The current function overwrite data in the dataframe.
        """
        data = []
        errors = 0
        for filename in filenames:
            df = self.parse_file(filename)
            if df is None:
                errors += 1
            else:
                data.append(df)

        temp_df = pd.concat(data, axis=1)
        temp_df = temp_df.transpose()
        print(temp_df.columns[0])
        if overwrite:
            self.__genes_counts = temp_df
            self.length = len(filenames) - errors
        else:
            self.__genes_counts = pd.concat([self.__genes_counts, temp_df], axis=0)
            self.length += (len(filenames) - errors)
    

        if self.__verbose:
            print(f"Files loaded. Action accomplished with {errors} inside {len(filenames)} files.")

    def load_annotations(self, file: str, overwrite = True) -> None:
        """
        Loads sample annotations from the XML file.

        Args:
            file (str): Path to the XML file containing sample annotations.
            overwrite (bool): The current function overwrite data in the dataframe.
        """
        xtree = et.parse(file)
        xroot = xtree.getroot()
        
        errors = 0

        temp_df = pd.DataFrame(columns = ["id","cns subregion", "tissue type", "sample group"])
        for child in xroot.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Sample"):
            try:
                temp_sample_id = child.attrib['iid']

                tags_map = {}
                for child2 in child.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Characteristics"):
                    if(child2.attrib["tag"] in self.__sample_annotations.columns):
                        tags_map[child2.attrib["tag"]] = child2.text.strip()
                temp_df = pd.concat([temp_df, pd.DataFrame({'id': [temp_sample_id], **tags_map})])
            except:
                errors += 1

        temp_df.set_index(keys='id', drop=True, inplace=True)
        if overwrite:
            self.__sample_annotations = temp_df
        else:
            self.__sample_annotations = pd.concat([self.__sample_annotations, temp_df], axis=0)
            
        if self.__verbose:
            print(f"Annotations loaded. Action accomplished with {errors} errors inside {len(xroot)} files.")

    def save_csv (self, filename: str) -> None:
        ''' 
        Save the dataframe __genes_counts to a CSV file.
        Args:
            filename (str): Path to the file where the dataframe will be saved.
        '''

        self.__genes_counts.to_csv(filename, index_label=self.__index_label)
        return

    def load_csv (self, filename: str) -> None:
        ''' 
        Load the dataframe __genes_counts from a CSV file.
        Args:
            filename (str): Path to the file where the dataframe will be loaded from.
        '''

        df = pd.read_csv(filename)
        df.rename(index=df[self.__index_label], inplace=True)
        df.drop(columns=self.__index_label, axis = 1, inplace= True)
        self.__genes_counts = df
        return
    
    #endregion
    
    #region Statistical Analysis
    
    def get_descriptive_statistics(self, by_gene: bool = True) -> pd.DataFrame:
        """
            If it's calculating desc statistics by gene, I do the logarithm transformation.
        """
        res = pd.DataFrame(columns = ["mean", "median", "std"])
        
        axis = 0
        if not by_gene:
            axis = 1

        res['mean'] = self.__genes_counts.mean(axis=axis)
        res['median'] = self.__genes_counts.median(axis=axis)
        res['std'] = self.__genes_counts.std(axis=axis)
        
        if by_gene:
            res['mean'] = np.log1p(res['mean'])
            res['median'] = np.log1p(res['median'])
            res['std'] = np.log1p(res['std'])

        return res

    def apply_log_transformation(self) -> pd.DataFrame:
        
        return

    #endregion

    #region Sample Description
    
    def get_sample_distribution_by_column(self, column: str) -> plt.Axes:
        return self.__sample_annotations[column].value_counts().plot(kind='bar', rot=45)
    
    #endregion
    
    #region RNA Counts
    
    #endregion

    def get_std_hist(self) -> plt.Axes:
        return self.__genes_counts.std(axis=1).plot.hist()
    
    def get_median_hist(self) -> plt.Axes:
        return self.__genes_counts.median(axis=1).plot.hist()
    
    def get_mean_hist(self) -> plt.Axes:
        return self.__genes_counts.mean(axis=1).plot.hist()
    
    def get_std(self) -> None:
        self.std_dev = self.__genes_counts.std()

        # Plot standard deviation
        plt.figure(figsize=(12, 6))

        plt.plot(self.std_dev, label='Standard Deviation')

        plt.legend(loc='best')
        plt.title('Standard Deviation of Each Gene')
        plt.xlabel('Gene')
        plt.ylabel('Value')

        plt.show()
        return
    
    def get_median(self) -> None:
        self.median = self.__genes_counts.median()

        # Plot median
        plt.figure(figsize=(12, 6))

        plt.plot(self.median, label='Median')

        plt.legend(loc='best')
        plt.title('Median of Each Gene')
        plt.xlabel('Gene')
        plt.ylabel('Value')

        plt.show()
        return
    
    def get_mean(self) -> None:
        # Plot mean
        self.mean = self.__genes_counts.mean()
        plt.figure(figsize=(12, 6))

        plt.plot(self.mean, label='Mean')

        plt.legend(loc='best')
        plt.title('Mean of Each Gene')
        plt.xlabel('Gene')
        plt.ylabel('Value')

        plt.show()
        return
        
    def get_mean_median_std(self) -> None:
        '''
        This method plots the mean, median and standard deviation of each gene and shows the plot.
        '''
        # mean = self.__genes_counts.mean()
        # median = self.__genes_counts.median()
        # std_dev = self.__genes_counts.std()

        # Plot mean, median, and standard deviation
        plt.figure(figsize=(12, 6))

        plt.plot(self.mean, label='Mean')
        plt.plot(self.median, label='Median')
        plt.plot(self.std_dev, label='Standard Deviation')

        plt.legend(loc='best')
        plt.title('Mean, Median, and Standard Deviation of Each Gene')
        plt.xlabel('Gene')
        plt.ylabel('Value')

        plt.show()
        return

    def log_transformation(self) -> pd.DataFrame:
        # Apply log transformation
        log_genes_counts = np.log1p(self.__genes_counts)

        mean = log_genes_counts.mean()
        median = log_genes_counts.median()
        std_dev = log_genes_counts.std()

        # Plot mean, median, and standard deviation
        plt.figure(figsize=(12, 6))

        plt.plot(mean, label='Mean')
        plt.plot(median, label='Median')
        plt.plot(std_dev, label='Standard Deviation')

        plt.legend(loc='best')
        plt.title('Mean, Median, and Standard Deviation of Each Gene (Log Transformed)')
        plt.xlabel('Gene')
        plt.ylabel('Value')

        plt.show()
        return log_genes_counts