set optimizer = 'sequential_pipe';

explain select lmin('5,10,1,4');
select lmin('5,10,1,4');
select lmin('5,10,2,4');

create table udf_lmin ( x string );
insert into udf_lmin values ('5,10,1,4');
insert into udf_lmin values ('2');
insert into udf_lmin values ('1,2,3');
insert into udf_lmin values ('10,2');
select * from udf_lmin;

explain select lmin(x) from udf_lmin;
select lmin(x) from udf_lmin;
