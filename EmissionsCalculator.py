
while(1):
    print "---------------Emissions Calculator-------------------------"
    desc = "This program calculates \
    emissions species concentratons from a combustion chamber using the CFD-CRN method\n\n"
    print desc
    print "Select one of the following options:\n"
    print "1) Read Fluent Casefile\n2) Read Datafile\n3) Generate Cluster\n4) Solve CRN\n5) Post Process"
    choice=input()
    if choice==1:
        import casefilepy
        casefilepy.convert("CFD_CRN.cas")
        continue
    if choice==2:
        import datfilepy
        datfilepy.readdata("CFD_CRN.dat")
        continue
    if choice==3:
        import CRN_Gen
        CRN_Gen.Generate()
        continue
    if choice==4:
        import ReactorGen_min
        ReactorGen_min.Gen(energy='off')
        continue
    if choice==5:
        import PostProc
        PostProc
        continue
    if choice==exit:
        break
    else:
        print "Invalid choice"

