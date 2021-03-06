/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *******************************************************************************/
package org.bgi.flexlab.gaea.data.structure.bam.filter;

import htsjdk.samtools.SAMRecord;
import org.bgi.flexlab.gaea.data.structure.bam.filter.util.FiltersMethod;
import org.bgi.flexlab.gaea.data.structure.bam.filter.util.MalformedReadFilter;
import org.bgi.flexlab.gaea.data.structure.bam.filter.util.ReadsFilter;
import org.bgi.flexlab.gaea.data.structure.bam.filter.util.SamRecordFilter;
import org.bgi.flexlab.gaea.data.structure.region.Region;

public class HaplotypeCallerFilter implements SamRecordFilter{
    HCReadsFilter readsFilter = new HCReadsFilter();

    MalformedReadFilter malformedReadFilter = new MalformedReadFilter();

    @Override
    public boolean filter(SAMRecord sam, Region region) {
        return (malformedReadFilter.filter(sam, region) || readsFilter.filter(sam, region) ||
                FiltersMethod.matchBaseNumberZero(sam) ||
                (region != null && !region.isSamRecordInRegion(sam)));
    }

    public class HCReadsFilter implements SamRecordFilter {

        @Override
        public boolean filter(SAMRecord sam, Region region) {
            return FiltersMethod.filterDuplicateRead(sam)
                    || FiltersMethod.filterMappingQualityUnavailable(sam,0)
                    || FiltersMethod.filterMappingQualityUnavailable(sam,255)
                    || FiltersMethod.filterUnmappedReads(sam)
                    || FiltersMethod.FailsVendorQualityCheckFilter(sam);
        }
    }

}
