import re
import numpy as np
import pandas as pd  # type: ignore
from sklearn.decomposition import PCA  # type: ignore
from sklearn.preprocessing import StandardScaler  # type: ignore
from sklearn.manifold import TSNE  # type: ignore
from sklearn.preprocessing import StandardScaler


class Samples:
    def __init__(self, verbose: bool = True, output_path: str = "outputs/genes_samples.csv") -> None:
        self.__genes_samples = pd.DataFrame()

        self.__stats = None
        self.__normalized_genes = pd.DataFrame()

        # PCA
        self.__reduced_samples = pd.DataFrame()
        self.__reduced_genes = pd.DataFrame()

        # tSNE - A set of (reduction, perplexity), in order to store multiple reductions with different perplexity.
        self.__reduced_genes_tSNE = []

        # configs
        self.__index_label = "sample"
        self.__output_path = output_path
        self.__verbose = verbose
        return

    @staticmethod
    def log_transform_df(df: pd.DataFrame) -> pd.DataFrame:
        return np.log1p(df)

    # region Imports

    def __save_csv(self) -> None:
        self.__genes_samples.to_csv(
            self.__output_path, index_label=self.__index_label)
        return

    def __load_csv(self) -> None:
        df = pd.read_csv(self.__output_path)
        df.rename(index=df[self.__index_label], inplace=True)
        df.drop(columns=self.__index_label, axis=1, inplace=True)
        self.__genes_samples = df
        return

    def __parse_file(self, filename: str) -> pd.DataFrame:
        try:
            df = pd.read_table(filename)
            sample_name = re.search("GSM\d+", filename).group()
            df.rename(index=df["gene/TE"], inplace=True)
            df.drop(columns=df.columns[0], axis=1, inplace=True)
            df.rename(columns={df.columns[0]: sample_name}, inplace=True)
            return df
        except:
            return None

    def load_data(self, filenames: list[str], check_from_csv=True) -> None:
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
                    print("Files loaded. Action accomplished with %d error(s) inside %d files.".format(
                        errors, len(filenames)))

        except:
            self.load_data(filenames, check_from_csv=False)
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

    # region PCA

    def __center_data(self, samples: pd.DataFrame) -> StandardScaler:
        scaler = StandardScaler()
        return scaler.fit_transform(samples)

    def __reduce_with_pca(self, df_to_reduce: pd.DataFrame, n_components=2) -> pd.DataFrame:
        pca = PCA(n_components=n_components)

        scaled_data = self.__center_data(df_to_reduce)
        pca = pca.fit(scaled_data)
        pca_data = pca.transform(scaled_data)

        vars = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
        labels = [f"PC{x}" for x in range(1, len(vars) + 1)]

        return pd.DataFrame(pca_data, index=df_to_reduce.index, columns=labels)

    def reduce_to_2d_per_gene(self, n_components=2):
        self.__reduced_genes = self.__reduce_with_pca(
            self.__genes_samples.T, n_components=n_components)
        return

    # endregion

    # region tSNE

    def __reduce_with_tSNE(self, df_to_reduce: pd.DataFrame, perplexity: int, n_components=2):
        reduced = TSNE(
            n_components=n_components,
            perplexity=perplexity
        ).fit_transform(df_to_reduce)

        return pd.DataFrame(reduced, index=df_to_reduce.index)

    def reduce_to_2d_per_gene_tSNE(self, perplexity: int, n_components=2):
        self.__reduced_genes_tSNE.append((self.__reduce_with_tSNE(
            self.__genes_samples.T, perplexity=perplexity, n_components=n_components), perplexity))
        return

    # endregion
