import xml.etree.ElementTree as et
import matplotlib.pyplot as plt
import pandas as pd
import re

# Here, I can create methods that give mean, std and median about the genes.
# I can also create a method called concat where, given more RNASeq instances
# I can return a new one (I create multiple instances of RNASeq (all the files)
# and then concatenate them together to show the difference???)
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
        self.__sample_annotations = pd.DataFrame(columns = ["id", "cns subregion", "tissue type", "sample group"])

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
        
    def load_single_data(self, filename: str) -> None:
        """
        Loads gene counts data from a single file.

        Args:
            filename (str): Path to the data file.
        """
        new_dataframe = self.parse_file(filename)
        new_dataframe = new_dataframe.transpose()
        self.__genes_counts = pd.concat([self.__genes_counts, new_dataframe], axis=1)
        self.length += 1
        
        if self.__verbose:
            if new_dataframe is None:
                print("Couldn't load the file. There have been errors.")
            else:
                print("File loaded.")
            
    def load_bulk_data(self, filenames: list[str]) -> None:
        """
        Loads gene counts data from multiple files.

        Args:
            filenames (list): List of paths to data files.
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
        self.__genes_counts = pd.concat([self.__genes_counts, temp_df], axis=1)
        self.length += (len(filenames) - errors)

        if self.__verbose:
            print(f"Files loaded. Action accomplished with {errors} inside {len(filenames)} files.")

    def load_annotations(self, file: str) -> None:
        """
        Loads sample annotations from the XML file.

        Args:
            file (str): Path to the XML file containing sample annotations.
        """
        xtree = et.parse(file)
        xroot = xtree.getroot()
        
        errors = 0

        for child in xroot.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Sample"):
            try:
                temp_sample_id = child.attrib['iid']

                tags_map = {}
                for child2 in child.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Characteristics"):
                    if(child2.attrib["tag"] in self.__sample_annotations.columns):
                        tags_map[child2.attrib["tag"]] = child2.text.strip()
                temp_df = pd.DataFrame({'id': [temp_sample_id], **tags_map})
                self.__sample_annotations = pd.concat([self.__sample_annotations, temp_df], axis=0)
            except:
                errors += 1
            
        if self.__verbose:
            print(f"Annotations loaded. Action accomplished with {errors} inside {len(xroot)} files.")

    #endregion
    
    #region Statistical Analysis

    def get_std_hist(self) -> plt.Axes:
        return self.__genes_counts.std(axis=1).plot.hist()
    
    def get_median_hist(self) -> plt.Axes:
        return self.__genes_counts.median(axis=1).plot.hist()
    
    def get_mean_hist(self) -> plt.Axes:
        return self.__genes_counts.mean(axis=1).plot.hist()
    
    def get_sample_distribution_by_column(self, column: str) -> plt.Axes:
        return self.__sample_annotations[column].value_counts().plot(kind='bar', rot=45)

    #endregion
    