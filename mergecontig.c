#include "mergecontig.h"
#include "string.h"

CtgDB* init_ctgdb(void ) {
	CtgDB *db;
	
	db = (CtgDB*)malloc(sizeof(CtgDB));
	db->ctgnum = 0;
	db->ctgs = init_ctglist(6);

	return db;
}

CtgDB* load_ctgdb(FileReader *fr1, FileReader *fr2) {
	uint32_t id = 0, i = 0;
	CtgDB *db;
	uuhash *map = init_uuhash(1023);
	uint32_t key, val;
	uuhash_t h;
	int len = 0;
	char *seq = NULL;
	String *line = init_string(1);

	while (fread_table(fr2) != -1) {
		key = atoi(get_col_str(fr2, 1));
		val = atoi(get_col_str(fr2, 4));
		h.key = key;
		h.val = val;

		if (!exists_uuhash(map, h)) {
			put_uuhash(map, h);
		}
	}

	db = (CtgDB*)malloc(sizeof(CtgDB));
	db->ctgnum = 0;
	db->ctgs = init_ctglist(6);

	while (fread_line(line, fr1) != -1) {
		if (line->string[0] == 'E') {
			if (len != 0) {
				Ctg contig;
				contig.id = id;
				contig.cls_id = i;
				h.key = id;
				h.val = 0;
				contig.old_clsid = get_uuhash(map, h)->val;
				contig.sz = 1;
				contig.seq = strdup(seq);
				db->ctgnum++;
				push_ctglist(db->ctgs, contig);
				i++;
			}
			free(seq); seq = NULL;
			len = 0;
			id = atoi(line->string+2);
		} else if (line->string[0] == 'S') {
			if (len < (int)strlen(line->string+2)) {
				len = (int)strlen(line->string+2);
				free(seq); seq = NULL;
				seq = strdup(line->string+2);
			}
		}
	}

	Ctg contig;
	contig.id = id;
	contig.cls_id = i;
	h.key = id;
	h.val = 0;
	contig.old_clsid = get_uuhash(map, h)->val;
	contig.sz = 1;
	contig.seq = strdup(seq);
	db->ctgnum++;
	push_ctglist(db->ctgs, contig);

	free(seq);
	free_string(line);
	free_uuhash(map);
	return db;
}

void print_ctgdb(CtgDB *db) {
	uint32_t i;
	Ctg *contig;

	for (i = 0; i < count_ctglist(db->ctgs); i++) {
		contig = ref_ctglist(db->ctgs, i);
		fprintf(stdout, "%d %d %d %d %s\n", contig->id, contig->cls_id, contig->old_clsid, contig->sz, contig->seq);
		fflush(stdout);
	}
	//fprintf(stdout, "%d\n", db->ctgnum);
}

void free_ctgdb(CtgDB *db) {
	//uint32_t i;
	//Ctg *contig;

	/*
	for (i = 0; i < count_ctglist(db->ctgs); i++) {
		contig = ref_ctglist(db->ctgs, i);
		free(contig->seq);
	}
	*/
	free_ctglist(db->ctgs);
	free(db);
}

void free_load_ctgdb(CtgDB *db) {
	uint32_t i;
	Ctg *contig;

	for (i = 0; i < count_ctglist(db->ctgs); i++) {
		contig = ref_ctglist(db->ctgs, i);
		free(contig->seq);
	}

	free_ctglist(db->ctgs);
	free(db);
}

int aln_cmp(const void *p0, const void *p1, void *ref) {
	PWcontig *t0, *t1;
	t0 = (PWcontig*)p0;
	t1 = (PWcontig*)p1;
	if (t0->score < t1->score)
		return 1;
	else if(t0->score > t1->score)
		return -1;
	else
		return 0;
	ref = ref;
}

PWDB* pw_aln_contigs(CtgDB *db) {
	uint32_t i, j, n;
	int k, mn, mm, off0, off1, aln_len;
	PWDB *pwdb;
	Ctg *c0, *c1;
	pwdb = (PWDB*)malloc(sizeof(PWDB));

	pwdb->pwctgs = init_pwctglist(6);
	pwdb->hp = init_heap(aln_cmp, pwdb);
	pwdb->ctgv = db->ctgs;
	AlnParam ap = {10, 2, 2, aln_sm_nt, 16, 75};

	n = db->ctgnum;
	
	for (i = 0; i < n-1; i++) {
		c0 = ref_ctglist(db->ctgs, i);
		for (j = i+1; j < n; j++) {
			c1 = ref_ctglist(db->ctgs, j);
			AlnAln *aa;
			mn = mm = 0;
			off0 = off1 = -1;
			aa = aln_stdaln(c0->seq, c1->seq, &ap, 0, 1);
			aln_len = strlen(aa->out1);
			for (k = 0; k < aln_len; k++) {
				if (aa->out1[k] == '-' || aa->out2[k] == '-') continue;
				if (aa->out1[k] != aa->out2[k]) mm++;
				mn++;
			}
			PWcontig *pwc = (PWcontig*)malloc(sizeof(PWcontig));
			pwc->id0 = c0->cls_id;
			pwc->id1 = c1->cls_id;
			pwc->overlap = mn;
			pwc->score = aa->score;
			pwc->het = (float)mm/mn; 
			push_heap(pwdb->hp, pwc);
			push_pwctglist(pwdb->pwctgs, pwc);
			//fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\t%.3f\n", c0->cls_id, c1->cls_id, pwc->id0, pwc->id1, mn, mm, pwc->het);
			//fprintf(stdout, "%s\n%s\n", c0->seq, c1->seq);
			//fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\n%s\n%s\n%s\n\n", aa->start1, aa->end1,aa->start2, aa->end2, pwc->score, pwc->overlap, aa->out1, aa->outm, aa->out2);

			//fprintf(stdout, "%s\n%s\n%s\n\n", aa->out1, aa->outm, aa->out2);
			fflush(stdout);
			aln_free_AlnAln(aa);
		}
	}

	return pwdb;
}
PWDB* clustering_ctg(PWDB *db, uint32_t min_overlap, float het) {
	PWDB *ret;
	ret = db;
	uint32_t i, j, p, q;
	
	PWcontig *poped;

	while ((poped = pop_heap(db->hp)) != NULL) {
		if (poped->overlap >= min_overlap && (poped->het - het <= 0)) {
		p = poped->id0;
		q = poped->id1;

		for (i = p; i != (ref_ctglist(ret->ctgv, i))->cls_id; i = (ref_ctglist(ret->ctgv, i))->cls_id)
			ref_ctglist(ret->ctgv, i)->cls_id = ref_ctglist(ret->ctgv, ref_ctglist(ret->ctgv, i)->cls_id)->cls_id;
		for (j = q; j != (ref_ctglist(ret->ctgv, j))->cls_id; j = (ref_ctglist(ret->ctgv, j))->cls_id)
			ref_ctglist(ret->ctgv, j)->cls_id = ref_ctglist(ret->ctgv, ref_ctglist(ret->ctgv, j)->cls_id)->cls_id;
		if (i == j) continue;
		if (ref_ctglist(ret->ctgv, i)->sz < ref_ctglist(ret->ctgv, j)->sz) {
			ref_ctglist(ret->ctgv, i)->cls_id = j;
			ref_ctglist(ret->ctgv, j)->sz += ref_ctglist(ret->ctgv, i)->sz;
		} else {
			ref_ctglist(ret->ctgv, j)->cls_id = i;
			ref_ctglist(ret->ctgv, i)->sz += ref_ctglist(ret->ctgv, j)->sz;
		}
		}
	}

	return ret;
}

int cmp_ctg_clsid(const void *p0, const void *p1) {
	Ctg *t0, *t1;
	t0 = (Ctg*)p0;
	t1 = (Ctg*)p1;
	if (t0->cls_id == t1->cls_id) return 0;
	if (t0->cls_id < t1->cls_id) return 1;
	return -1;
}

void print_clusters(PWDB *db) {
	int last_cid = -1;
	int line_num = 0;
	uint32_t i;
	Ctg *ctg;

	ctglist *t = db->ctgv;
	
	qsort(as_array_ctglist(t), count_ctglist(t), sizeof(Ctg), cmp_ctg_clsid);
	for (i = 0; i < count_ctglist(t); i++) {
		ctg = ref_ctglist(t, i);
		if (last_cid != (int)ctg->cls_id && line_num > 0)
			printf("\n");
		last_cid = ctg->cls_id;
		line_num++;
		printf("%d ", ctg->id);
	}
	printf("\n");

	return;

}

void execute_pwaln(CtgDB *db, uint32_t min_overlap, float het, uint32_t max_nctg) {
	PWDB *pwaln;
	CtgDB *tdb = init_ctgdb();

	uint32_t i;
	int last_oldcid = -1;
	int line_num = 0;
	Ctg *ctg;
	int id = 0;

	
	for (i = 0; i < count_ctglist(db->ctgs); i++) {
		ctg = ref_ctglist(db->ctgs, i);
		if (last_oldcid != (int)ctg->old_clsid) {
			if (line_num > 0) {
				tdb->ctgnum = count_ctglist(tdb->ctgs);
				if (tdb->ctgnum > 1 && tdb->ctgnum <= max_nctg) { //magic number 50
					//print_ctgdb(tdb);
					
					pwaln = pw_aln_contigs(tdb);
					pwaln = clustering_ctg(pwaln, min_overlap, het);
					print_clusters(pwaln);
					free_pwdb(pwaln); 
					free_ctgdb(tdb);
				} else {
					free_ctgdb(tdb);
				}
				id = 0;
				tdb = init_ctgdb();
			}
			last_oldcid = (int)ctg->old_clsid;
			line_num++;
			ctg->cls_id = id;
			push_ctglist(tdb->ctgs, *ctg);
			id++;
		} else {
			line_num++;
			ctg->cls_id = id;
			push_ctglist(tdb->ctgs, *ctg);
			id++;
		}
	}
	tdb->ctgnum = count_ctglist(tdb->ctgs);
	if (tdb->ctgnum > 1 && tdb->ctgnum <= max_nctg) { //magic number 50
		//print_ctgdb(tdb);
		
		pwaln = pw_aln_contigs(tdb);
		pwaln = clustering_ctg(pwaln, min_overlap, het);
		print_clusters(pwaln);
		free_pwdb(pwaln);
		free_ctgdb(tdb);
	} else {
		free_ctgdb(tdb);
	}
}

void free_pwdb(PWDB *db) {
	uint32_t i;
	PWcontig *pw;

	for (i = 0; i < count_pwctglist(db->pwctgs); i++) {
		pw = get_pwctglist(db->pwctgs, i);
		free(pw);
	}
	free_pwctglist(db->pwctgs);
	free_heap(db->hp);
	
	free(db);
}
