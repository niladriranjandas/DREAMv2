
def genAtomSite(selection="known",
                model=None,
                occupancies=None,
                bfactors=None,
                coordPrecision=3,
                addICode=False,
                dataBlock=None,
                ):
    """
    Given the <m atomSel>.AtomSel selection, return a CifDatablock
    with an atom_site category added.

    Fields not populated are simply not written.

    If specified, occupancies and bfactors must be arrays with the same length
    as selection.

    The number of digits after the decimal of the output coordinates is given
    by coordPrecision.

    If dataBlock is not specified, a new CifDatablock is created.

    """

    if dataBlock==None:
        from cif import CifDatablock
        dataBlock = CifDatablock()
        pass
    from cif import CifCategory
    cat = CifCategory()
    
    from selectTools import convertToAtomSel
    selection=convertToAtomSel(selection)

    #
    atomTypes = sorted( list(set([atom.atomName()[0] for atom in selection])) )
    cat.addKey("symbol")
    for type in atomTypes:
        cat.addValue("symbol",type)
        pass
    dataBlock["atom_type"] = cat

    cat = CifCategory()
    for key in ("group_PDB",         #          ATOM
                "id",                #                 1-offset index
                "type_symbol",       #      required
                "label_atom_id",     #      atomName
                "label_alt_id",      #      altLoc
                "label_comp_id",     #      residueName
                "label_asym_id",     #      segmentName
                "label_entity_id",   #    . (or 1)
                "label_seq_id",      #       residueNum
                "Cartn_x",           #      specify digits after decimal point 
                "Cartn_y", 
                "Cartn_z", 
                "auth_comp_id", # .
                "auth_asym_id", # .
                "auth_seq_id", # .
                ):
        cat.addKey(key)
        pass
    if addICode:    cat.addkey("pdbx_PDB_ins_code") #    iCode
    if occupancies: cat.addKey("occupancy")
    if bfactors:    cat.addKey("B_iso_or_equiv")
    if model:       cat.addKey("pdbx_PDB_model_num")


    for i,atom in enumerate(selection):
        pos=atom.pos()
        cat.addValue("group_PDB", "ATOM")
        cat.addValue("id", str(i+1))              
        cat.addValue("type_symbol",     atom.atomName()[0])
        cat.addValue("label_atom_id",   atom.atomName())
        cat.addValue("label_alt_id",   ".")
        cat.addValue("label_comp_id",   atom.residueName())
        cat.addValue("label_asym_id",   atom.segmentName())
        cat.addValue("label_entity_id", ".")
        cat.addValue("label_seq_id",    str(atom.residueNum()))
        if addICode:
            cat.addValue("pdbx_PDB_ins_code", ".")
            pass
        cat.addValue("Cartn_x",         "%.*f" % (coordPrecision,pos[0]))
        cat.addValue("Cartn_y",         "%.*f" % (coordPrecision,pos[1]))
        cat.addValue("Cartn_z",         "%.*f" % (coordPrecision,pos[2]))
        if occupancies:
            cat.addValue("occupancy",       "%.*f" % (occupancyPrecision,
                                                      occupancies[i])    )
            pass
        if bfactors:
            cat.addValue("B_iso_or_equiv",  "%.*f" % (bfactorPrecision,
                                                      bfactors[i])     )
            pass
        if addICode:
            cat.addValue("pdbx_PDB_ins_code", ".")
            pass
#        cat.addValue("auth_atom_id", # . 
        cat.addValue("auth_comp_id",   atom.residueName())
        cat.addValue("auth_asym_id",   atom.segmentName())
        cat.addValue("auth_seq_id",   str(atom.residueNum()))
        if model:
            cat.addValue("pdbx_PDB_model_num", str(model) )
            pass
        pass
    dataBlock["atom_site"] = cat
    return dataBlock
    
                         














        


#loop_
#_atom_site.group_PDB          ATOM
#_atom_site.id                 1-offset index
#_atom_site.type_symbol        .
#_atom_site.label_atom_id      atomName
#_atom_site.label_alt_id       .
#_atom_site.label_comp_id      residueName
#_atom_site.label_asym_id      segmentName
#_atom_site.label_entity_id    . (or 1)
#_atom_site.label_seq_id       residueNum
#_atom_site.pdbx_PDB_ins_code  .
#_atom_site.Cartn_x            specify digits after decimal point 
#_atom_site.Cartn_y 
#_atom_site.Cartn_z 
#_atom_site.occupancy 
#_atom_site.B_iso_or_equiv 
#_atom_site.Cartn_x_esd 
#_atom_site.Cartn_y_esd 
#_atom_site.Cartn_z_esd 
#_atom_site.occupancy_esd 
#_atom_site.B_iso_or_equiv_esd 
#_atom_site.pdbx_formal_charge 
#_atom_site.auth_seq_id 
#_atom_site.auth_comp_id 
#_atom_site.auth_asym_id 
#_atom_site.auth_atom_id 
#_atom_site.pdbx_PDB_model_num 
#ATOM 1     N N    . MET A 1 1  ? -14.136 1.321   3.616   1.00 0.93 ? ? ? ? ? ? 1  MET A N    1  
