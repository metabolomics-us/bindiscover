class VennTableQuery():
    def build_query(
        self,
        dropdown_triplet_selection_value,
        slider_percent_present_value,
        toggle_average_true_value,
        radio_items_filter_value
    ):
    
        sod_parallel_list=[element.split(' - ') for element in dropdown_triplet_selection_value]
        if toggle_average_true_value==True:
            intensity_type='intensity_average'
        elif toggle_average_true_value==False:
            intensity_type='intensity_median'
        slider_percent_present_value*=0.01
    
        main_line_list=[
            f'case when (max(percent_present) filter (where "species"=\'{sod_parallel_list[i][0]}\' and "organ"=\'{sod_parallel_list[i][1]}\' and "disease"=\'{sod_parallel_list[i][2]}\')>{slider_percent_present_value}) then (max({intensity_type}) filter (where "species"=\'{sod_parallel_list[i][0]}\' and "organ"=\'{sod_parallel_list[i][1]}\' and "disease"=\'{sod_parallel_list[i][2]}\')) else null end as "{dropdown_triplet_selection_value[i]}",\n' for i in range(len(dropdown_triplet_selection_value))
        ]

        main_line=''.join(main_line_list)
        main_line=main_line[:-2]+'\n'

        if radio_items_filter_value=='no_filter':
            radio_filter_string=''
        elif radio_items_filter_value=='common':
            radio_filter_string_list=[
                f'(\"{element}\" is not null) ' for element in dropdown_triplet_selection_value
            ]
            radio_filter_string=' and '.join(radio_filter_string_list)
            radio_filter_string=' where '+radio_filter_string

        self.query=f'''
            select * from (
            select * from (
            select
            bin,
            compound,\n'''+main_line+'''from
            non_ratio_table nrt
            group by
            bin,compound 
            order by
            bin) as foo'''+radio_filter_string+'''
			) as bar inner join 
			non_ratio_table_fix nrtf on
			bar.bin=nrtf.compound_identifier            
            '''
