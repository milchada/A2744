(for ((i=0; i<10; i++)); do
echo >&2 "$i"
convert +append s1-s2_Data_00000$i*png s1-s2_Data_00000$i-Slice_density.png
done)
