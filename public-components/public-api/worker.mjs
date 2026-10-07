import compounds from '../public-data/dataset.json' with {type:'json'};
import foods from '../public-data/foods-pilot.json' with {type:'json'};
import components from '../COMPONENTS.json' with {type:'json'};
import provenance from './provenance.json' with {type:'json'};
import openapi from './openapi.json' with {type:'json'};
import {createHandler} from './handler.mjs';
export default {fetch:createHandler({compounds,foods,components,provenance,openapi})};
