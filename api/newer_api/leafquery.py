class LeafQuery():
    
    def __init__(self,
        triplet_from,
        triplet_to,
        bin_type
    ):

        if bin_type=='known':
            self.query=f'''
                select compound_id,fold_change_average,significance_welch from differential_analysis_knowns dak where
                (dak.triplet_from={triplet_from}) and (dak.triplet_to={triplet_to})
                '''
        elif bin_type=='unknown':
            self.query=f'''
                select compound_id,fold_change_average,significance_welch from differential_analysis_unknowns dau where
                (dau.triplet_from={triplet_from}) and (dau.triplet_to={triplet_to})
                '''
        elif bin_type=='class':
            self.query=f'''
                select compound_id,fold_change_average,significance_welch from differential_analysis_classes dac where
                (dac.triplet_from={triplet_from}) and (dac.triplet_to={triplet_to})
                '''
        #print(self.query)