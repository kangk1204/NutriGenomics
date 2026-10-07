"""Read-only food composition with confirmed portions and original units."""
import math
from pathlib import Path


class FoodRegistry:
    def __init__(self,root,database):
        self.root=Path(root).resolve();self.database=Path(database).resolve()
        if not self.database.is_relative_to(self.root):raise ValueError('Food database must remain inside research root')

    def connection(self):
        if not self.database.is_file():raise FileNotFoundError('Frozen food database is not installed')
        import duckdb
        connection=duckdb.connect(str(self.database),read_only=True)
        connection.execute('SET threads=2')
        return connection

    def search(self,query='',limit=25):
        if not isinstance(query,str) or len(query)>200:raise ValueError('Food query must be at most 200 characters')
        limit=max(1,min(int(limit),100))
        with self.connection() as c:
            rows=c.execute('SELECT food_id,source_id,source_food_id,description,data_type,publication_date FROM food WHERE contains(lower(description),lower(?)) ORDER BY food_id LIMIT ?',[query,limit]).fetchall()
        columns=['food_id','source_id','source_food_id','description','data_type','publication_date']
        return {'items':[dict(zip(columns,r)) for r in rows],'unit':'unique_food_record','query':query}

    def composition(self,food_id,grams):
        if not math.isfinite(grams) or not 0<grams<=10000:raise ValueError('Confirmed grams must be positive and at most 10000')
        with self.connection() as c:
            food=c.execute('SELECT food_id,source_id,source_food_id,description,data_type,publication_date FROM food WHERE food_id=?',[food_id]).fetchone()
            if food is None:raise KeyError('Unknown food record')
            rows=c.execute('SELECT nutrient_id,nutrient_name,amount,unit,basis,source_row_id,source_id FROM food_composition WHERE food_id=? ORDER BY nutrient_id,source_row_id',[food_id]).fetchall()
        valid,excluded=[],0
        for nutrient,name,amount,unit,basis,row_id,source in rows:
            if basis!='source100g' or amount is None or not math.isfinite(float(amount)) or float(amount)<0 or not unit:
                excluded+=1;continue
            # Source exports contain impossible mass values. Preserve them in the
            # ETL rejection ledger and never serve more than 100 g mass per 100 g food.
            factor={'G':1.,'MG':.001,'UG':.000001,'MCG':.000001}.get(unit.upper())
            if factor is not None and float(amount)*factor>100:
                excluded+=1;continue
            valid.append({'nutrient_id':nutrient,'nutrient_name':name,'source_amount':float(amount),'source_basis':basis,
                          'portion_amount':float(amount)*grams/100.,'unit':unit,'source_row_id':row_id,'source_id':source})
        return {'food':dict(zip(['food_id','source_id','source_food_id','description','data_type','publication_date'],food)),
                'confirmed_grams':grams,'nutrients':valid,'excluded_incompatible_rows':excluded,
                'calculation':'source amount per 100 g × confirmed grams / 100',
                'clinical_effect_inferred':False,'image_estimated_portion':False}
