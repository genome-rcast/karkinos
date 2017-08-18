/*
Copyright Hiroki Ueda

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
package jp.ac.utokyo.rcast.karkinos.summary;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

public class ChromBand {

	Connection memorycon = null;
	String tablename = "cbdata";

	public ChromBand(String csv) {

		loadCsv(csv);

	}

	public ChromBand(String csv, String tbname) {

		tablename = tbname;
		loadCsv(csv);

	}

	public void loadCsv(String csv) {

		System.out.println("start loadng");

		try {

			Class.forName("org.h2.Driver");
			memorycon = DriverManager.getConnection("jdbc:h2:mem:mydb", "sa",
					"");
			String sqlcreate = "create table " + tablename
					+ " as select * from csvread('" + csv + "',null,null,'\t')";
			PreparedStatement crs = memorycon.prepareStatement(sqlcreate);
			crs.execute();
			System.out.println("end loadng csv");

		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	List<CNAInterval> list = null;

	public List<CNAInterval> getList() throws SQLException {

		List<CNAInterval> list = new ArrayList<CNAInterval>();

		String sql = "select * from " + tablename;
		PreparedStatement crs = memorycon.prepareStatement(sql);

		ResultSet rs = crs.executeQuery();
		while (rs.next()) {

			CNAInterval ci = new CNAInterval();
			ci.chr = rs.getString(1);
			ci.start = rs.getInt(2);
			ci.end = rs.getInt(3);
			ci.name = rs.getString(4);
			ci.setGieStain(rs.getString(5));
			list.add(ci);

		}
		this.list = list;
		return list;

	}

	public String getBand(String chr, int start, int end) throws SQLException {

		StringBuffer sb = new StringBuffer();
		//
		if (list == null) {
			getList();
		}
		for (CNAInterval ci : list) {
			if (!ci.chr.equals(chr)) {
				continue;
			}
			if (ci.intercect(start, end)) {
				if (sb.length() > 0) {
					sb.append(",");
				}
				sb.append(ci.name);
			}
		}

		return sb.toString();

	}

	public void close() throws SQLException {
		if (memorycon != null) {
			try {
				memorycon.createStatement().execute("drop table " + tablename);
			} catch (Exception ex) {
			}
			;
			if (!memorycon.isClosed()) {
				memorycon.close();
			}
		}
		memorycon = null;
	}
}
