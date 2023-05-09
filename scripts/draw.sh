#!/bin/bash

useLog=0
cal="calib/ADCstep.root"
dr="recon.: all peaks within 8 ns";
tw="-1400 9800 32"
suf="coin8ns.defaultOffset"

trig90D="trig.: > 90 mV in T0L/R"; trig90="trigAB90mV90mV"
trig9070D="trig.: > 90 mV in T0L/R, < 70 mV in T1/2"; trig9070="trigABbCbD90mV90mV70mV70mV"
trig4f90D="trig.: > 90 mV in T0L/R, > 25 mV in T1/2"; trig4f90="trigABCD90mV90mV25mV25mV"
trig4f40D="trig.: > 40 mV in T0L/R, > 25 mV in T1/2"; trig4f40="trigABCD40mV40mV25mV25mV"
trig5f35D="trig.: > 35 mV in T0L/R, > 25 mV in T1/2"; trig5f35="trigABCD35mV35mV25mV25mV"
#trig30D="trig.: > 30 mV in T0L/R"; trig30="trigAB30mV30mV"
trig25D="trig.: > 25 mV in T0L/R"; trig25="trigAB25mV25mV"
trig20D="trig.: > 20 mV in T0L/R"; trig20="trigAB20mV20mV"
trig30D="trig.: > 30 mV in T0L/R"; trig30="trigT0_30mV"
trig40D="trig.: > 30 mV in T0L/R"; trig40="trigT0_40mV"

dtype="coin" # 4-fold, 2-fold, 3-fold, single

#for ct in "foundTrigger&&trig_n==4" "foundTrigger&&trig_n==3&&trig_hD==0" "foundTrigger&&trig_n==3&&trig_hC==0" "foundTrigger&&trig_n==2"
for ct in "foundTrigger&&trig_n0==2&&trig_hC==0&&trig_hD==0" "foundTrigger&&trig_n0==2&&(trig_hC>0||trig_hD>0)" "foundTrigger"
do
    if [ ${ct} == "foundTrigger" ]
    then
        dt="all types of triggers";
        nt="trigall"
    elif [ ${ct} == "foundTrigger&&trig_n0==2&&(trig_hC>0||trig_hD>0)" ]
    then
        dt="triggers w/ peaks in T0 and T1/T2";
        nt="trigT0T12"
    elif [ ${ct} == "foundTrigger&&trig_n0==2&&trig_hC==0&&trig_hD==0" ]
    then
        dt="triggers w/ 2 peaks in T0, no peak in T1/T2";
        nt="trigT0"
    fi
#    for ce in "evt_n==4" "evt_n==3&&evt_hD==0" "evt_n==3&&evt_hC==0" "evt_n==3&&evt_n0==1" "evt_n==2&&evt_n0==0" "evt_n==2&&evt_hC==0&&evt_hD==0" "evt_n==1&&evt_hC>0" "evt_n==1&&evt_hD>0" "evt_n==1&&evt_n0>0"
    for ce in "evt_n0==0&&evt_hC>0"
    do
        if [ ${ce} == "evt_n==4" ]
        then
            de="events with 4-fold coincidence"
            ne="evt4-fold"
        elif [ ${ce} == "evt_n0==0&&evt_hC>0" ]
        then
            de="events with T1 hits, no hits in T0"
            ne="evtT1T12"
        elif [ ${ce} == "evt_n==3&&evt_hD==0" ]
        then
            de="events with 3-fold coincidence, excl. T2"
            ne="evt3-foldmD"
        elif [ ${ce} == "evt_n==3&&evt_hC==0" ]
        then
            de="events with 3-fold coincidence, excl. T1"
            ne="evt3-foldmC"
        elif [ ${ce} == "evt_n==3&&evt_n0==1" ]
        then
            de="events with 3-fold coincidence, excl. T0L/R"
            ne="evt3-foldm0"
        elif [ ${ce} == "evt_n==2&&evt_n0==0" ]
        then
            de="events with 2-fold coincidence in T12"
            ne="evt2-fold12"
        elif [ ${ce} == "evt_n==2&&evt_hC==0&&evt_hD==0" ]
        then
            de="events with 2-fold coincidence in T0"
            ne="evt2-fold0"
        elif [ ${ce} == "evt_n==1&&evt_hC>0" ]
        then
            de="events with single hit in T1"
            ne="evtsingle1"
        elif [ ${ce} == "evt_n==1&&evt_hD>0" ]
        then
            de="events with single hit in T2"
            ne="evtsingle2"
        elif [ ${ce} == "evt_n==1&&evt_n0>0" ]
        then
            de="events with single hit in T0L/R"
            ne="evtsingle0"
        fi
        ce="${ce}&&evt_t!=trig_t&&(evt_h0==0||evt_h0>30)&&(evt_hC==0||evt_hC>30)&&(evt_hD==0||evt_hD>30)" # TODO: this is temporary. To be tuned
#        for sizeChar in 'w' 'h' 'a'
        for sizeChar in 'h'
        do
            pre="thrA20mVB20mVC10mVD10mV.evtall"
#            run="run00370-01000";./BinaryFiles/bin/draw -D -F -r "${run}.${pre}${nt}.${ne}" -d "${run}, ${dr}, ${trig9070D}" -s ${sizeChar} -c "${dt}" -o "${de}" -S "${ct} ${ce}" -E ${dtype} -C ${cal} -T "${tw}" -L ${useLog} anaout/run00370-00505.${pre}${trig9070}.event.root anaout/run00584-00624.${pre}${trig9070}.event.root anaout/run00773-01000.${pre}${trig9070}.event.root
#            run="run00509-00530";./BinaryFiles/bin/draw -D -F -r "${run}.${pre}${nt}.${ne}" -d "${run}, ${dr}, ${trig4f40D}" -s ${sizeChar} -c "${dt}" -o "${de}" -S "${ct} ${ce}" -E ${dtype} -C ${cal} -T "${tw}" -L ${useLog} anaout/run00509-00517.${pre}${trig4f40}.event.root anaout/run00518-00530.${pre}${trig4f40}.event.root
#            run="run00531-00583";./BinaryFiles/bin/draw -D -F -r "${run}.${pre}${nt}.${ne}" -d "${run}, ${dr}, ${trig5f35D}" -s ${sizeChar} -c "${dt}" -o "${de}" -S "${ct} ${ce}" -E ${dtype} -C ${cal} -T "${tw}" -L ${useLog} anaout/run00531-00583.${pre}${trig5f35}.event.root
#            run="run00700-00772";./BinaryFiles/bin/draw -D -F -r "${run}.${pre}${nt}.${ne}" -d "${run}, ${dr}, ${trig9070D}" -s ${sizeChar} -c "${dt}" -o "${de}" -S "${ct} ${ce}" -E ${dtype} -C ${cal} -T "${tw}" -L ${useLog} anaout/run00700-00772.${pre}${trig9070}.event.root
#            run="run01001-01040";./BinaryFiles/bin/draw -D -F -r "${run}.${pre}${nt}.${ne}" -d "${run}, ${dr}, ${trig9070D}" -s ${sizeChar} -c "${dt}" -o "${de}" -S "${ct} ${ce}" -E ${dtype} -C ${cal} -T "${tw}" -L ${useLog} anaout/run01001-01006.${pre}${trig9070}.event.root anaout/run01007-01040.${pre}${trig9070}.event.root
#            run="run01436-01507";./BinaryFiles/bin/draw -D -F -r "${run}.${pre}${nt}.${ne}" -d "${run}, ${dr}, ${trig25D}"   -s ${sizeChar} -c "${dt}" -o "${de}" -S "${ct} ${ce}" -E ${dtype} -C ${cal} -T "${tw}" -L ${useLog} anaout/run01436-01507.${pre}${trig25}.event.root
#            run="run01508-01529";./BinaryFiles/bin/draw -D -F -r "${run}.${pre}${nt}.${ne}" -d "${run}, ${dr}, ${trig25D}"   -s ${sizeChar} -c "${dt}" -o "${de}" -S "${ct} ${ce}" -E ${dtype} -C ${cal} -T "${tw}" -L ${useLog} anaout/run01508-01529.${pre}${trig25}.event.root
            pre="thrA20mVB20mVC20mVD20mVE20mVF20mV.evtall"
            TW="1300 1400 0.8"
#            run="run01530-01599";./BinaryFiles/bin/draw -D -F -r "${run}.${pre}${nt}.${ne}" -d "${run}, ${dr}, ${trig25D}"   -s ${sizeChar} -c "${dt}" -o "${de}" -S "${ct} ${ce}" -E ${dtype} -C ${cal} -T "${tw}" -L ${useLog} anaout/${run}.${pre}${trig25}.event.root
#            run="run01600-01699";./BinaryFiles/bin/draw -D -F -r "${run}.${pre}${nt}.${ne}" -d "${run}, ${dr}, ${trig25D}"   -s ${sizeChar} -c "${dt}" -o "${de}" -S "${ct} ${ce}" -E ${dtype} -C ${cal} -T "${tw}" -L ${useLog} anaout/${run}.${pre}${trig25}.event.root
#            run="run01700-01798";./BinaryFiles/bin/draw -D -F -r "${run}.${pre}${nt}.${ne}" -d "${run}, ${dr}, ${trig25D}"   -s ${sizeChar} -c "${dt}" -o "${de}" -S "${ct} ${ce}" -E ${dtype} -C ${cal} -T "${tw}" -L ${useLog} anaout/${run}.${pre}${trig25}.event.root
#            run="run01802-01881";./BinaryFiles/bin/draw -D -F -r "${run}.${pre}${nt}.${ne}" -d "${run}, ${dr}, ${trig25D}"   -s ${sizeChar} -c "${dt}" -o "${de}" -S "${ct} ${ce}" -E ${dtype} -C ${cal} -T "${tw}" -L ${useLog} anaout/${run}.${pre}${trig25}.event.root
#            run="run01530-01881";./BinaryFiles/bin/draw -D -F -r "${run}.${pre}${nt}.${ne}" -d "${run}, ${dr}, ${trig25D}"   -s ${sizeChar} -c "${dt}" -o "${de}" -S "${ct} ${ce}" -E ${dtype} -C ${cal} -T "${tw}" -L ${useLog} anaout/run01600-01699.${pre}${trig25}.event.root anaout/run01530-01599.${pre}${trig25}.event.root anaout/run01700-01798.${pre}${trig25}.event.root anaout/run01802-01881.${pre}${trig25}.event.root
#             run="run02327-02462";echo ./BinaryFiles/bin/draw -D -F \"60 8600\" -r \"${run}.${pre}.${trig40}.${suf}.${nt}.${ne}\" -d \"${run}, ${dr}, ${trig40D}\"   -s ${sizeChar} -c \"${dt}\" -o \"${de}\" -S \"${ct} ${ce}\" -E ${dtype} -C ${cal} -t \"${tw}\" -T \"${TW}\" -L ${useLog} anaout/${run}.${pre}.${trig40}.${suf}.event.root 
#             run="run01530-01882";echo ./BinaryFiles/bin/draw -D -F \"60 8600\" -r \"${run}.${pre}.${trig40}.${suf}.${nt}.${ne}\" -d \"${run}, ${dr}, ${trig40D}\"   -s ${sizeChar} -c \"${dt}\" -o \"${de}\" -S \"${ct} ${ce}\" -E ${dtype} -C ${cal} -t \"${tw}\" -T \"${TW}\" -L ${useLog} anaout/${run}.${pre}.${trig40}.${suf}.event.root 
             run="run02720-02893";echo ./BinaryFiles/bin/draw -D -F \"60 8600\" -r \"${run}.${pre}.${trig40}.${suf}.${nt}.${ne}\" -d \"${run}, ${dr}, ${trig40D}\"   -s ${sizeChar} -c \"${dt}\" -o \"${de}\" -S \"${ct} ${ce}\" -E ${dtype} -C ${cal} -t \"${tw}\" -T \"${TW}\" -L ${useLog} anaout/${run}.${pre}.${trig40}.${suf}.event.root 
        done
    done
done
