for i in "$@"
do
case $i in
	-r)
	RIBO=true
	;;
	-s)
	STRUCT=true
	;;
	-c)
	CM=true
	;;
	--train)
	TYPE='train'
	;;
	--holdout)
	TYPE='holdout'
	;;
	--test)
	TYPE='test'
	;;
esac
done

if [ $STRUCT ]; then
readlink -f ../input_data/StructureData/${TYPE}/*.bpseq
fi

if [ $RIBO ]; then
readlink -f ../input_data/RiboswitchData/${TYPE}/*.bpseq
fi

if [ $CM ]; then
readlink -f ../input_data/ChemMappingData/${TYPE}/*.bpseq
fi
