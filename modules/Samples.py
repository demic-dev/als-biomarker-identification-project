import re
import numpy as np
import pandas as pd  # type: ignore
from sklearn.decomposition import PCA # type: ignore
from sklearn.preprocessing import StandardScaler # type: ignore
from sklearn.manifold import TSNE # type: ignore
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler



class Samples:
    def __init__(self, verbose: bool = True, output_path: str = "outputs/genes_samples.csv") -> None:
        self.__genes_samples = pd.DataFrame()
        
        self.__stats = None

        self.__reduced_samples = pd.DataFrame()
        self.__reduced_genes = pd.DataFrame()
        
        self.__reduced_samples_tSNE = pd.DataFrame()
        self.__reduced_genes_tSNE = pd.DataFrame()
        
        self.__index_label = "sample"
        self.__output_path = output_path
        self.__verbose = verbose
        self.__normalized_genes = pd.DataFrame()
        return

    # region Imports

    def __save_csv (self) -> None:
        self.__genes_samples.to_csv(self.__output_path, index_label=self.__index_label)
        return

    def __load_csv (self) -> None:
        df = pd.read_csv(self.__output_path)
        df.rename(index=df[self.__index_label], inplace=True)
        df.drop(columns=self.__index_label, axis = 1, inplace= True)
        self.__genes_samples = df
        return

    def __parse_file(self, filename: str) -> pd.DataFrame:
        try:
            df = pd.read_table(filename)
            sample_name = re.search("GSM\d+", filename).group()
            df.rename(index= df["gene/TE"], inplace=True)
            df.drop(columns=df.columns[0], axis=1, inplace=True)
            df.rename(columns={ df.columns[0]: sample_name }, inplace = True)
            return df
        except:
            return None
        
    def load_data(self, filenames: list[str], check_from_csv = True) -> None:
        try:
            if check_from_csv:
                self.__load_csv()
            else:
                data = []
                errors = 0
                for filename in filenames:
                    df = self.__parse_file(filename)
                    if df is None:
                        errors += 1
                    else:
                        data.append(df)

                temp_df = pd.concat(data, axis=1)
                temp_df = temp_df.transpose()

                self.__genes_samples = temp_df
                self.__save_csv()

                if self.__verbose:
                    print(f"Files loaded. Action accomplished with {errors} error(s) inside {len(filenames)} files.")
            
        except:
            self.load_data(filenames, check_from_csv = False)
            pass

    # endregion

    # region Setters and Getters 
    
    def get_samples(self) -> pd.DataFrame:
        return self.__genes_samples
    
    def get_stats(self):
        return self.__stats
    
    def get_reduced_samples(self) -> pd.DataFrame:
        return self.__reduced_samples
    
    def get_reduced_genes(self) -> pd.DataFrame:
        return self.__reduced_genes
    
    def get_reduced_samples_tSNE(self) -> pd.DataFrame:
        return self.__reduced_samples_tSNE

    def get_reduced_genes_tSNE(self) -> pd.DataFrame:
        return self.__reduced_genes_tSNE

    def toggle_verbose(self) -> None:
        self.__verbose = not self.__verbose
        return

    def is_verbose(self) -> bool:
        return self.__verbose
    
    def get_normalized_genes(self) -> pd.DataFrame:
        return self.__normalized_genes

    # endregion
    
    # region Statistical Analysis

    def get_basic_statistics(self) -> pd.DataFrame:
        self.__stats = self.__genes_samples.describe()
        return self.__stats 
    
    def get_mean_median_std_dev(self) -> pd.DataFrame:
        mean = self.__stats.loc['mean']
        median = self.__stats.loc['50%']  # Median is the 50th percentile
        std_dev = self.__stats.loc['std']

        # Create a new figure
        plt.figure(figsize=(15, 10))

        # Subplot for mean
        plt.subplot(2, 2, 1)
        plt.plot(mean)
        plt.title('Mean')

        # Subplot for median
        plt.subplot(2, 2, 2)
        plt.plot(median)
        plt.title('Median')

        # Subplot for standard deviation
        plt.subplot(2, 2, 3)
        plt.plot(std_dev)
        plt.title('Standard Deviation')

        # Subplot for all
        plt.subplot(2, 2, 4)
        plt.plot(mean, label='Mean')
        plt.plot(median, label='Median')
        plt.plot(std_dev, label='Standard Deviation')
        plt.title('All Statistics')
        plt.legend()

        # Display the plots
        plt.tight_layout()
        plt.show()
        
        return

    def get_descriptive_statistics(self, by_gene: bool = True) -> pd.DataFrame:
        columns = ["mean", "median", "std"]
        res = pd.DataFrame(columns=columns)
        
        axis = 1
        if by_gene:
            axis = 0

        res[columns[0]] = self.__genes_samples.mean(axis=axis)
        res[columns[1]] = self.__genes_samples.median(axis=axis)
        res[columns[2]] = self.__genes_samples.std(axis=axis)
        
        if by_gene:
            res[columns[0]] = np.log1p(res[columns[0]])
            res[columns[1]] = np.log1p(res[columns[1]])
            res[columns[2]] = np.log1p(res[columns[2]])

        return res
    
    def get_normalized_genes_data(self) -> pd.DataFrame:

        # Assuming df is your DataFrame
        scaler = StandardScaler()
        self.__normalized_genes = pd.DataFrame(scaler.fit_transform(self.__genes_samples), columns=self.__genes_samples.columns)
        return self.__normalized_genes
    # endregion
    
    # region PCA
    
    def __center_data(self, samples: pd.DataFrame) -> StandardScaler:
        scaler = StandardScaler()
        return scaler.fit_transform(samples)
    
    def __reduce_to_2d(self, df_to_reduce: pd.DataFrame) -> pd.DataFrame:
        pca = PCA(n_components= 2)

        scaled_data = self.__center_data(df_to_reduce)
        pca = pca.fit(scaled_data)
        pca_data = pca.transform(scaled_data)

        vars = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
        labels = [f"PC{x}" for x in range(1, len(vars) + 1)]

        return pd.DataFrame(pca_data, index=df_to_reduce.index, columns = labels)
    
    def reduce_to_2d_per_sample(self):
        self.__reduced_samples = self.__reduce_to_2d(self.__genes_samples)
        return
    
    def reduce_to_2d_per_gene(self):
        self.__reduced_genes =self.__reduce_to_2d(self.__genes_samples.T)
        return
    
    # endregion
    
    # region tSNE
    
    def __reduce_to_2d_with_tSNE(self, df_to_reduce: pd.DataFrame):
        scaled_data = self.__center_data(df_to_reduce)
        reduced = TSNE(
            n_components=2,
            learning_rate='auto',
            init='random',
            perplexity=20
        ).fit_transform(scaled_data)

        return pd.DataFrame(reduced, index = df_to_reduce.index)
        
    def reduce_to_2d_per_sample_tSNE(self):
        self.__reduced_samples_tSNE = self.__reduce_to_2d_with_tSNE(self.__genes_samples)
        return
    
    def reduce_to_2d_per_gene_tSNE(self):
        self.__reduced_genes_tSNE = self.__reduce_to_2d_with_tSNE(self.__genes_samples.T)
        return
    
    # endregion

