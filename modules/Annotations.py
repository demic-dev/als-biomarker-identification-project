import xml.etree.ElementTree as et
import pandas as pd  # type: ignore


class Annotations:
    def __init__(self, verbose: bool = True, output_path="outputs/annotations.csv") -> None:
        self.__columns = ["cns subregion", "sample group"]
        self.__sample_annotations = pd.DataFrame(self.__columns)

        self.__index_label = "sample"
        self.__output_path = output_path
        self.__verbose = verbose
        return

    def __save_csv(self) -> None:
        self.__sample_annotations.to_csv(self.__output_path)
        return

    def __load_csv(self, filename: str) -> None:
        df = pd.read_csv(filename)
        df.rename(index=df[self.__index_label], inplace=True)
        df.drop(columns=self.__index_label, axis=1, inplace=True)
        self.__sample_annotations = df
        return

    def load_annotations(self, file: str, check_from_csv=True) -> None:
        try:
            if check_from_csv:
                self.__load_csv()
            else:
                xtree = et.parse(file)
                xroot = xtree.getroot()

                errors = 0

                temp_df = pd.DataFrame(columns=self.__columns)
                for child in xroot.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Sample"):
                    try:
                        temp_sample_id = child.attrib['iid']

                        tags_map = {}
                        for child2 in child.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Characteristics"):
                            if (child2.attrib["tag"] in self.__columns):
                                tags_map[child2.attrib["tag"]
                                         ] = child2.text.strip()
                        temp_df = pd.concat(
                            [temp_df, pd.DataFrame(tags_map, index=[temp_sample_id])])
                    except:
                        errors += 1

                # temp_df.set_index(keys='id', drop=True, inplace=True)

                self.__sample_annotations = temp_df
                self.__save_csv()

                if self.__verbose:
                    print("Annotations loaded. Action accomplished with %d errors inside %d files.".format(errors, len(xroot)))
        except Exception as e:
            self.load_annotations(file, check_from_csv=False)
            pass

    # region Setters and Getters

    def get_annotations(self) -> pd.DataFrame:
        return self.__sample_annotations

    def get_annotations_column(self, column: str) -> pd.DataFrame:
        if column not in self.__columns:
            authorised_values = ", ".join(
                [f"`{column}`" for column in self.__columns])
            raise ValueError(
                "The column argument should have one of these values: %s".format(authorised_values)
            )
        return self.__sample_annotations[column]

    # endregion

    # def get_descriptive_statistics(self, by_gene: bool = True) -> pd.DataFrame:
    #     res = pd.DataFrame(columns = ["mean", "median", "std"])

    #     axis = 0
    #     if not by_gene:
    #         axis = 1

    #     res['mean'] = self.__genes_counts.mean(axis=axis)
    #     res['median'] = self.__genes_counts.median(axis=axis)
    #     res['std'] = self.__genes_counts.std(axis=axis)

    #     if by_gene:
    #         res['mean'] = np.log1p(res['mean'])
    #         res['median'] = np.log1p(res['median'])
    #         res['std'] = np.log1p(res['std'])

    #     return res
