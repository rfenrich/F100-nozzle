
for f in *.pdf
do
		pdfcrop $f
		ff="${f%.*}-crop.pdf"
		echo $ff
		mv $ff $f
done


