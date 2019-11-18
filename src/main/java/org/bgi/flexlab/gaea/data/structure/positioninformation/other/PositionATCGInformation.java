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
package org.bgi.flexlab.gaea.data.structure.positioninformation.other;

import org.bgi.flexlab.gaea.data.structure.bam.SAMInformationBasic;
import org.bgi.flexlab.gaea.data.structure.positioninformation.CalculatePositionInforamtionInterface;
import org.bgi.flexlab.gaea.data.structure.positioninformation.CompoundInformation;
import org.bgi.flexlab.gaea.data.structure.positioninformation.PositionInformationUtils;

public class PositionATCGInformation extends PositionInformationUtils<ATCGCount> implements CalculatePositionInforamtionInterface<SAMInformationBasic>{
	
	public PositionATCGInformation(int windowSize) {
		super(windowSize);
	}
	
	@SuppressWarnings("rawtypes")
	@Override
	public void add(CompoundInformation posInfo) {
		// TODO Auto-generated method stub
		int base = posInfo.getBinaryBase();
		if(base < 4) {
			info.get(posInfo.distBetweenRefPosAndWinStart()).add(1, base);
		}
	}
	
	
}
