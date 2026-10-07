from nutriomics_evidence.graph import connect,node,structure_coverage

def test_available_structure_does_not_certify_ctd_cross_database_identity(tmp_path):
    db=connect(tmp_path/'graph.sqlite')
    node(db,'InChIKey:XLYOFNOQVPJJNP-UHFFFAOYSA-N','chemical','water',{'full_inchikey':'XLYOFNOQVPJJNP-UHFFFAOYSA-N'})
    node(db,'MESH:D014867','chemical','Water')
    with db:
        coverage=structure_coverage(db)
    assert coverage=={'chemical_nodes':2,'chemical_nodes_with_full_structure':1,
                      'chemical_nodes_without_full_structure':1,'cross_database_exact_structure_mappings':0}
    db.close()
