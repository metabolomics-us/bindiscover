class SunburstQuery():

    def init():
        pass

    def build_query(
        self,
        compound
    ):

        self.query=f'''
            select 
            species,
            organ,
            disease,
            max(intensity_average) filter (where "bin"={compound}) as intensity_average,
            max(intensity_median) filter (where "bin"={compound}) as intensity_median,
            max(percent_present) filter (where "bin"={compound}) as percent_present 
            from 
            non_ratio_table nrt 
            group by
            species, organ, disease
            order by 
            species,organ,disease  
        '''
