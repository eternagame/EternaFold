./generate_datalist.sh -s --holdout > s_holdout.txt
./generate_datalist.sh -c --holdout > c_holdout.txt
./generate_datalist.sh -r --holdout > r_holdout.txt
./generate_datalist.sh -s -c --holdout > sc_holdout.txt
./generate_datalist.sh -s -r --holdout > sr_holdout.txt
./generate_datalist.sh -c -r --holdout > cr_holdout.txt
./generate_datalist.sh -s -c -r --holdout > scr_holdout.txt
./generate_datalist.sh -s --test > s_test.txt
./generate_datalist.sh -c --test > c_test.txt
./generate_datalist.sh -r --test > r_test.txt
./generate_datalist.sh -s -c --test > sc_test.txt
./generate_datalist.sh -s -r --test > sr_test.txt
./generate_datalist.sh -c -r --test > cr_test.txt
./generate_datalist.sh -s -c -r --test > scr_test.txt
./generate_datalist.sh -s --train > s_train.txt
./generate_datalist.sh -c --train > c_train.txt
./generate_datalist.sh -r --train > r_train.txt
./generate_datalist.sh -s -c --train > sc_train.txt
./generate_datalist.sh -s -r --train > sr_train.txt
./generate_datalist.sh -c -r --train > cr_train.txt
./generate_datalist.sh -s -c -r --train > scr_train.txt
