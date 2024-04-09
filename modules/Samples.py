import re
import numpy as np
import pandas as pd  # type: ignore

class Samples:
    def __init__(self, verbose: bool = True, output_path: str = "outputs/genes_samples.csv") -> None:
        self.__genes_samples = pd.DataFrame()
        
        self.__index_label = "sample"
        self.__output_path = output_path
        self.__verbose = verbose
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

    def toggle_verbose(self) -> None:
        self.__verbose = not self.__verbose
        return

    def is_verbose(self) -> bool:
        return self.__verbose

    # endregion
    
    # region Statistical Analysis

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
    
    # endregion