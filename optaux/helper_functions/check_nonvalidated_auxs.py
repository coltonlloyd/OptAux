import cobra

met_to_rs = {'EX_pydam_e': ['PDX5PS', 'PYDXK', 'PYDXNK'],
             'EX_orot_e': ['DHORTS', 'UPPRT', 'URIK2'],
             'EX_thr__L_e': ['PTHRpp', 'THRS'],
             'EX_pro__L_e': ['AMPTASEPG', 'P5CR'],
             'EX_skm_e': ['DHQTi'],
             'EX_cys__L_e': ['AMPTASECG', 'CYSS']}

for m, rs in met_to_rs.items():
    ijo = cobra.io.load_json_model(
        '/home/sbrg-cjlloyd/Desktop/ecoli_M_models/iML1515.json')
    ijo.reactions.EX_o2_e.lower_bound = -20

    ijo.reactions.BIOMASS_Ec_iML1515_core_75p37M.lower_bound = .1
    ijo.reactions.BIOMASS_Ec_iML1515_core_75p37M.upper_bound = .1

    for r in rs:
        for g in [i.id for i in ijo.reactions.get_by_id(r).genes]:
            print(ijo.genes.get_by_id(g).name,
                  [i.id for i in ijo.genes.get_by_id(g).reactions])
            ijo.genes.get_by_id(g).remove_from_model()

    ijo.objective = m
    ijo.reactions.get_by_id(m).lower_bound = -10
    ijo.optimize()
    print(m, ijo.solution.f)
    ijo.reactions.get_by_id(m).lower_bound = 0